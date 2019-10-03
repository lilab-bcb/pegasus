=============================================
How to run a containerized VM on Google Cloud
=============================================



Create a VM instance with a docker image
========================================

1. Go to `Google Cloud Platform console`_. Click **Navigation menu** at the top left corner. Click **Compute Engine** under *COMPUTE*. 

2. Click **VM instances** from the left and then click **CREATE INSTANCE** at the top. 

3. Choose **New VM instance** from the left and then fill out the form.

4. In the form, fill the **Name** field with the name of your choice, such as ``container``. Under **Machine type**, choose the number of CPUs and RAMs you want to use , e.g. *64 vCPUs*.

5. Click *Deploy a container image to this VM instance* under **Container**. Then fill in the docker image. For images deposited to Docker Hub, the image name is the same as what you will use in ``docker pull`` command. For example, fill in ``sccloud/sccloud:0.9.1`` to install pegasus image. 

6. Click **Advanced container options**. Change *Restart policy* to *On failure* and check *Allocate a buffer for STDIN* and *Allocate a pseudo-TTY*.

7. To mount host directory to docker container, click *Add volume* under **Volume mounts** and fill in the table: *Volume Type* = Directory, *Mount path* = /home/userid, *Host path* = /home/userid, *Mode* = Read/write. Then click *Done*. Note that if the host path does not exist, docker will not be installed correctly.

8. Click the *Change* buttion under **Boot disk** if you need to increase the disk space. The default disk space is only *10G*. For example, if you want to increase the boot disk space to *200*GB, change the number under *Size (GB)* to *200* and click *Select* at the bottom right corner.

9. Click **Create** at the bottom of the form.



Run the created VM
==================

1. You should be able to see the VM name with a green check marker under **VM instances**. Click *SSH* under the *Connect* column to launch a command-line shell. 

2. Check if the docker image is successfully installed. If so, you should be able to see a welcome message (in organce color) telling you that you need to use ``docker attach`` to access your containers. If you cannot find the welcome message, it is possible that the cloud engine is still working on deploying the docker image. Wait for one minute and proceed to the next step.

3. Type ``docker ps`` several times until you see your docker image appears (e.g. sccloud/sccloud:0.9.1 under the IMAGE column). Copy the corresponding container name under the NAMES column (e.g. klt-instance-3-uwtu).

4. Attach to the container using the following command::

	docker attach container_name

5. Use *ctrl+p, ctrl+q* to detach from the container and use ``docker attach container_name`` to reattach to the container.

6. You can check the CPU info outside the container using the command below::

	cat /proc/cpuinfo



SSH remotely
============

1. Add your public key to the VM instance by following the *Adding or removing instance-level public SSH keys* section in `Managing SSH Keys in Metadata`_.

2. Use the command below to ssh remotely to the VM instance (also see *Connecting using third-party tools* section in `Connecting to Instances Using Advanced Methods`_)::
	
	ssh -i [PATH_TO_PRIVATE_KEY] [USERNAME]@[EXTERNAL_IP_ADDRESS]


.. _Google Cloud Platform console: https://console.cloud.google.com
.. _Managing SSH Keys in Metadata: https://cloud.google.com/compute/docs/instances/adding-removing-ssh-keys#edit-ssh-metadata
.. _Connecting to Instances Using Advanced Methods: https://cloud.google.com/compute/docs/instances/connecting-advanced#provide-key
