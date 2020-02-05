.. _Python (3.5+): https://www.python.org/downloads/
.. _Perl (5.26+): https://www.perl.org/get.html
.. _Blast (2.2+): https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _GitHub repository: https://github.com/SternLabTAU/AccuNGS/

Downloading and installing AccuNGS
==================================
AccuNGS code is present on our `GitHub repository`_. 
The different folders correspond to the different stages of AccuNGS as described in :doc:`overview`. 
AccuNGS has a few external dependencies:

* `Perl (5.26+)`_
* `Blast (2.2+)`_
* `Python (3.5+)`_

After downloading the code from the github repository, it is required to add the directory of the
installation to the ``base_calling/base_call_and_freqs.pl`` script (line 17).