.. _Pandas: https://pandas.pydata.org/
.. _Numpy: https://numpy.org/
.. _Scipy: https://www.scipy.org/

Variant calling
===============
This part is taking advantage of a homogeneous control population sequenced in parallel to the samples of interest. 
This homogeneous control acts as a marker for the baseline error levels of the protocol, 
which can be modeled by several gamma distributions that correspond to different erroneous substitutions (i.e. A>G, G>A, etc.)

Prerequisits
^^^^^^^^^^^^
This python script requires `Pandas`_, `Numpy`_ and `Scipy`_ python packages installed. 

Input
^^^^^

The inputs for the variant calling step are a FREQS file of the sample of interest, and the FREQS file of the homogenous control sample. 
They don't have to share the same genetic background (loci or organism), but the should share the same library preperation, 
sequencing and output generation methods.

===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
sample_freqs          Text           FREQS file of sample of interest
--------------------- -------------- --------------------------------
control_freqs         Text           FREQS file of a homogeneous control
--------------------- -------------- --------------------------------
-c / --coverage       Integer        Minimal sequencing depth to consider positions. Defaults to 1,000
--------------------- -------------- --------------------------------
-o / --output         Text           Output variant file. Defaults to "output.var.csv"
===================== ============== ================================

Invokation:

.. code-block:: bash

  python variant_calling/variant_caller.py ${sample_freqs} ${control_freqs} ${min_coverage} ${output_file}

Output
^^^^^^

The output by this step is similar to the output by 
`V-phaser2 <https://www.broadinstitute.org/viral-genomics/v-phaser-2>`_. 
Here's an :download:`example variants file <examples/example_output.var.csv>`.
Each possible variant for each loci with a high-enough coverage is assigned 
a line in the output file. The ``pval`` column can be used to filter variants, 
to allow a specific false-positive rate in variant calling. 
1% filtering is a reasonable filtering value. 

Example
^^^^^^^
For the sake of this example, we run on the :download:`example FREQS output <examples/example_output.freqs>`
on the :doc:`basecall` step on itself, using it as both a sample and a control. 
If the sample is not too diverse, it may be reasonable since extreme outliers (i.e. real variants 
and not sequencing errors) are not used to fit the gamma distributions.

.. code-block:: bash

  python variant_caller/variant_caller.py merged.fasta.freqs merged.fasta.freqs -c 500 -o ${output_file}

The result of this line of code is similar to this :download:`example variants file <examples/example_output.var.csv>`
