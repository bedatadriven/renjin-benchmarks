
Clinical eslII (Essential Statistics Learning 2nd Edition)
==========================================================

In clinica studies, information about disease state and progression of a large
number of patients is collected together with patient information that could
(technically) be used for early diagnosis and or prognosis of disease. In this
workflow clinical prostate cancer data of 97 patients is used to identify
molecular markers that correlate with disease stage. The analysis is mainly
inspired by the book '`Essential Statistical Learning`_'.

Analysis used in this workflow are cross-validation for Lasso penalized
regression fit and best predictive variable identification ('ncvreg' and
'leaps' packages); linear regression and fitting with L1 constraints ('stats'
and 'lasso2' packages); Lasso penalized Least Angle Regression with cross
validation ('lars' package); and fitting General Linear Model ('stats'
package).


Packages and Dependencies
-------------------------

There are 6 packages used in this workflow, which depend
on 2 additional packages from CRAN (dependencies)

**Used packages:**

* *CRAN*: ncvreg, boot, lars, lasso2, mda, leaps

**Package dependencies:**

* *CRAN*: class, MASS

Data
-------

source datasets from http://cran.r-project.org/web/packages/lasso2/lasso2.pdf

* Prostate dataset:
    * Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer.
    * Stamey, T., et al. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. Journal of Urology, 16: 1076-1083.
    
.. _Essential Statistical Learning: http://statweb.stanford.edu/~tibs/ElemStatLearn/

License
-------

* Copyright (c) 2015 Ieuan Clay based on code from `genbench <https://github.com/biolion/genbench>`_
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html

