.. _Pandas: https://pandas.pydata.org/
.. _Numpy: https://numpy.org/
.. _Scipy: https://www.scipy.org/

Variant calling
===============
This stage takes advantage of a homogeneous control population sequenced in parallel to the samples of interest. 
AccuNGS uses the homogeneous control sample to learn how process errors look like and fit gamma distribution to each type of mutation (e.g. A>G, G>A, etc.). 
It then assigns a p-value for all possible observed mutations obtained from the base-calling stage, 
representing the probability that the observed mutation frequency belongs to the corresponding gamma distribution.

Prerequisites
^^^^^^^^^^^^^
This python script requires `Pandas`_, `Numpy`_ and `Scipy`_ python packages installed. 

Input
^^^^^

The inputs for the variant calling stage are a FREQS file of the sample of 
interest, and the FREQS file of the homogenous control sample (see 
:ref:`basecall_output`). Ideally, they share the same genetic background (loci 
or organism), but most importantly they should share the same library 
preparation, sequencing and output generation methods.

===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
sample_freqs          Text           FREQS file of sample of interest
--------------------- -------------- --------------------------------
control_freqs         Text           FREQS file of a homogeneous control
--------------------- -------------- --------------------------------
``-c / --coverage``   Integer        Minimal sequencing depth to consider positions. Defaults to 1,000
--------------------- -------------- --------------------------------
``-o / --output``     Text           Output variant file. Defaults to ``output.var.csv``
===================== ============== ================================

Invokation:

.. code-block:: bash

  python variant_calling/variant_caller.py ${sample_freqs} ${control_freqs} ${min_coverage} ${output_file}

Output
^^^^^^

The output of this step is similar to the output by 
`V-phaser2 <https://www.broadinstitute.org/viral-genomics/v-phaser-2>`_. 
Here's an :download:`example variants file <examples/example_output.var.csv>`.
Each possible variant for each locus with a high-enough coverage is assigned 
a line in the output file. The ``pval`` column can be used to filter variants, 
to allow a specific false-positive rate in variant calling. 
During AccuNGS calibration, we've used p-value<=1% filtering value. 

.. note::
    Currently insertions or deletions are not supported in this variants file. 
    Use the output of the :doc:`basecall` step to identify such variants.

Example
^^^^^^^
For the sake of this example, we run on the :download:`example FREQS output <examples/example_output.freqs>`
from the :doc:`basecall` step using it as both a sample and a control for itself. 
If the sample is not too diverse, this actually might be reasonable since extreme outliers (i.e. real variants 
and not sequencing errors) are not used to fit the gamma distributions.

.. code-block:: bash

  python variant_caller/variant_caller.py merged.fasta.freqs merged.fasta.freqs -c 500 -o ${output_file}

The result of this line of code is similar to this :download:`example variants file <examples/example_output.var.csv>`.
