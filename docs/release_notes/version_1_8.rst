1.8.1 :small:`August 23, 2023`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Bug fix in cell marker JSON files for ``infer_cell_types`` function.

1.8.0 :small:`July 21, 2023`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Feature and Improvement**

* Updata ``human_immune`` and ``human_lung`` marker sets.
* Add ``mouse_liver`` marker set.
* Add `split_one_cluster <./api/pegasus.split_one_cluster.html>`_ function to subcluster one cluster into a specified number of subclusters.
* Update **neighbors** function to set ``use_cache=False`` by default, and adjust K to ``min(K, int(sqrt(n_samples)))``. [PR `272 <https://github.com/lilab-bcb/pegasus/pull/272>`_]
* In **infer_doublets** function, argument ``manual_correction`` now accepts a float number threshold specified by users for cut-off. [PR `275 <https://github.com/lilab-bcb/pegasus/pull/275>`_]

**Bug Fix**

* Fix divide by zero issue in ``integrative_nmf`` function. [PR `258 <https://github.com/lilab-bcb/pegasus/pull/258>`_]
* Compatibility with Pandas v2.0. [PR `261 <https://github.com/lilab-bcb/pegasus/pull/261>`_]
* Allow ``infer_doublets`` to use any count matrix with key name specified by users. [PR `268 <https://github.com/lilab-bcb/pegasus/pull/268>`_ Thanks to `Donghoon Lee <https://github.com/hoondy>`_]
