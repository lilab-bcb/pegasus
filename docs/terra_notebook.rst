Use ``pegasus`` on Terra Notebook
----------------------------------

You need to first have a `Terra <https://app.terra.bio/>`_ account.

1. Start Notebook Runtime on Terra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first time when you use Terra notebook, you need to create a runtime environment. 

On top-right panel of your workspace, click the following button within red circle:

.. image:: images/create_runtime.png
   :scale: 70 %
   :align: center

In the pop-out dialog (see image below), select ``CUSTOM ENVIRONMENT`` tab. Then in ``Container image`` field, choose one of Docker images from `Here <https://hub.docker.com/repository/docker/cumulusprod/pegasus-terra>`_. All the tags are for different versions of ``Pegasus``. For example, if you want to use Pegasus version ``0.16.1``, you should type ``cumulusprod/pegasus-terra:0.16.1``:

.. image:: images/runtime_setting.png
   :scale: 50 %
   :align: center

After that, set the computing resources you want to use in the new runtime environment within ``COMPUTE POWER`` field. And click ``REPLACE`` button, then ``CREATE`` button.

After waiting for 1-2 minutes, your runtime environment should be ready. In the pop-up dialog, click ``APPLY`` button to start.

This runtime environment is associated with your workspace, not just the current workspace. Therefore, after creation, you can start the same environment in your workspace by clicking the following button within red circle on top-right panel:


2. Create Your Terra Notebook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In ``NOTEBOOKS`` tab of your workspace, you can either create a blank notebook, or upload your local notebook. After creation, click ``EDIT`` button to enter the edit mode, and Terra will automatically start your notebook runtime.

When the start-up is done, you can type the following code to check if ``pegasus`` can be loaded and if it's the correct version you want to use::

	import pegasus as pg
	pg.__version__


3. Load Data into Runtime
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use your data on cloud (i.e. in the Google Bucket of your workspace), you should copy it into your notebook runtime with Google Cloud SDK::

	!gsutil -m cp gs://link-to-h5-file .

where ``gs://link-to-h5-file`` is the Google Bucket URL to your count matrix data file. 

After that, you can use Pegasus function to load it into memory.

4. Stop Notebook Runtime
^^^^^^^^^^^^^^^^^^^^^^^^^

When you are done with interactive analysis, don't forget to stop your notebook runtime by clicking the following button of the top-right panel of your workspace within red circle:

Otherwise, Google Cloud will keep on charging you. 