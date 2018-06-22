from firecloud import api as fapi
import argparse
import json
import subprocess
import os
import warnings

warnings.filterwarnings('ignore', 'Your application has authenticated', UserWarning, 'google')


def firecloud_run(method, config_name, workspace, wdl_inputs, method_version, upload):
    if type(wdl_inputs) != dict:
        with open(args.wdl_inputs, 'r') as f:
            inputs_json = json.loads(f.read())
            # Filter out any key/values that contain #, and escape strings with quotes as MCs need this to not be treated as expressions
            inputs = {k: "\"{}\"".format(v) for k, v in inputs_json.items() if '#' not in k}

    else:
        inputs = wdl_inputs
    sep = workspace.index('/')
    workspace_namespace = workspace[0:sep]
    workspace_name = workspace[sep + 1:]

    sep = method.index('/')
    method_namespace = method[0:sep]
    method_name = method[sep + 1:]

    config_namespace = workspace_namespace
    if config_name is None:
        config_name = method_name + '_config'

    ws = fapi.get_workspace(workspace_namespace, workspace_name)
    if ws.status_code == 404:
        ws = fapi.create_workspace(workspace_namespace, workspace_name)
        if ws.status_code != 201:
            raise ValueError('Unable to create workspace')
        workspace_id = ws.json()['workspaceId']
    else:
        workspace_id = ws.json()['workspace']['workspaceId']

    if upload is not None:
        for path in upload:
            path = os.path.abspath(path)
            if not os.path.exists(path):
                raise ValueError(path + ' does not exist.')

            subprocess.check_call(['gsutil', '-m', 'rsync', '-r', f, 'gs://' + workspace_id])

    if method_version is None:
        version = -1
        methods = fapi.list_repository_methods(method_name).json()
        for method in methods:
            if method['namespace'] == method_namespace:
                version = max(version, method['snapshotId'])
        if version == -1:
            raise ValueError(method_name + ' not found')

        method_version = version

    root_entity = None
    launch_entity = None
    method_body = {
        'name': config_name,
        'namespace': config_namespace,
        'methodRepoMethod': {'methodNamespace': method_namespace, 'methodName': method_name,
                             'methodVersion': method_version},
        'rootEntityType': root_entity,
        'prerequisites': {},
        'methodConfigVersion': 1,
        'deleted': False,
        'inputs': inputs,
        'outputs': {}
    }
    config_exists = fapi.get_workspace_config(workspace_namespace, workspace_name, config_namespace, config_name)

    if config_exists.status_code == 200:
        config_submission = fapi.update_workspace_config(workspace_namespace, workspace_name, config_namespace,
                                                         config_name, method_body)
        if config_submission.status_code != 200:
            raise ValueError('Unable to update workspace config')

    else:
        config_submission = fapi.create_workspace_config(workspace_namespace, workspace_name, method_body)
        if config_submission.status_code != 201:
            raise ValueError('Unable to create workspace config')

    launch_submission = fapi.create_submission(workspace_namespace, workspace_name, config_namespace, config_name,
                                               launch_entity, root_entity, "")
    if launch_submission.status_code == 201:
        submission_id = launch_submission.json()['submissionId']
        print('https://portal.firecloud.org/#workspaces/{}/{}/monitor/{}'.format(workspace_namespace, workspace_name,
                                                                                 submission_id))
    else:
        print(launch_submission.json())
        raise ValueError('Unable to launch submission - ' + launch_submission.json())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Upload files to a Google bucket and run a WDL')
    parser.add_argument('-m', '--method', dest='method', action='store', required=True,
                        help='Method name (e.g. regevlab/cellranger)')
    parser.add_argument('-c', '--config_name', dest='config_name', action='store', required=False,
                        help='Method configuration name')
    parser.add_argument('-w', '--workspace', dest='workspace', action='store', required=True,
                        help='Workspace name (foo/bar). The workspace will be created if it does not exist')
    parser.add_argument('-i', '--wdl-inputs', dest='wdl_inputs', action='store', required=True, help='WDL input JSON')
    parser.add_argument('-v', '--method_version', dest='method_version', action='store',
                        help='Method version. If not specified the latest version of the method will be used.')
    parser.add_argument('-u', '--upload', dest='upload', action='append',
                        help='Files or directories to upload to the workspace. Directories will be uploaded recursively.')
    args = parser.parse_args()
    firecloud_run(args.method, args.config_name, args.workspace, args.wdl_inputs, args.method_version, args.upload)
