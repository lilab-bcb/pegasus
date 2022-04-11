1.6.0 :small:`April 10, 2022`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Features**

* Add support for `scVI-tools <https://github.com/scverse/scvi-tools>`_:

    * Function `pegasus.run_scvi <./api/pegasus.run_scvi.html>`_, which is a wrapper of scVI for data integration.
    * Add a dedicated section for scVI method in Pegasus `batch correction tutorial <https://pegasus-tutorials.readthedocs.io/en/latest/_static/tutorials/batch_correction.html>`_.
    * Function `pegasus.train_scarches_scanvi <./api/pegasus.train_scarches_scanvi.html>`_ and `pegasus.predict_scarches_scanvi <./api/pegasus.predict_scarches_scanvi.html>`_ to wrap `scArches <https://github.com/theislab/scarches>`_ for transfer learning on single-cell data labels.
