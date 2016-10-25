
Clinical eslII (Essential Statistics Learning 2nd Edition)
==========================================================

In clinical studies, information about disease state and progression of a large
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

.. graphviz::
   :caption: Workflow to identify prostate cancer markers using clinical data.

   digraph CLINICAL_eslII {
		{rank=same pros_regsubset pros_l1ce pros_lars pros_cvlars pros_glm}
		DataP [shape = invhouse, label = "data(prostate)"];
		DataPTr [shape = invhouse, label = "Training set"];
		DataPTe [shape = invhouse, label = "Test set"];

		vs_ncvreg [shape = box; label = "Coordinate descent     \ncv.ncvreg()"];
		pros_regsubset  [shape = box; label = "regsubset()  "];
		pros_l1ce [shape = box; label = "l1ce()   "];
		pros_lars [shape = box; label = "lars()   "];
		pros_cvlars [shape = box; label = "cv.lars()   "];
		pros_glm [shape = box; label = "glm()    "];
		pros_cvglm [shape = box; label = "cv.glm()"];
		pros_rbind [label = "merge results    \ndo.call(rbind)"];
		pros_print [label = "print()"];

		DataP -> DataPTr;
		DataP -> DataPTe;
        DataP -> vs_ncvreg -> pros_rbind -> pros_print;
        DataPTr -> pros_regsubset -> pros_rbind;
        DataPTr -> pros_l1ce -> pros_rbind;
        DataPTr -> pros_lars -> pros_rbind;
        DataPTr -> pros_cvlars -> pros_rbind;
        DataPTr -> pros_glm -> pros_rbind;
        pros_glm -> pros_cvglm -> pros_rbind;
	}



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

