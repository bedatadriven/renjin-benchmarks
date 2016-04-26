********************
Functions
********************

This workflows contains many functions for plotting, identification of input data and computations.
You can read about each function input, what it does, and the resulting output.

++++++++++++++++++++
Plotting Functions
++++++++++++++++++++

coverAndKeyPlot
+++++++++++++++++
Creates a cover- and index-page with used naming/grouping.

Usage
^^^^^^^
coverAndKeyPlot(description, refName, WIDTH, HEIGHT)

Arguments
^^^^^^^^^^

description
     ... (Default is NULL).

refName
     Used in final tables/plots file name (default is "").
WIDTH
    defines the page width (default is 1000).
HEIGHT
    defines the page height (default is 1414).

Return
^^^^^^^
Initialization of a cover page with the given specifictions.

QCtablePlot
++++++++++++++++
Plotting based on simpleaffy QC function.

Usage
^^^^^^
QCtablePlot(Data, quality, sprep, lys, samplePrep, ratio, hybrid,
            percPres, bgPlot, scaleFact, WIDTH, HEIGHT, POINTSIZE)

Arguments
^^^^^^^^^^
Data
    ...
quality
    ... Default is NULL.
sprep
    ... default is NULL.
lys
    ... default is NULL.
samplePrep
    ... default is TRUE.
ratio
    ... default is TRUE.
hybrid
    ... default is TRUE.
percPres
    ... default is TRUE.
bgPlot
    ...
scaleFact
    ...
WIDTH
    ...
HEIGHT
    ...
POINTSIZE
    ...

Return
^^^^^^^
.. image:: QCtable.png

samplePrepPlot
++++++++++++++++++
Plots samplePrep controls

Usage
^^^^^^
samplePrepPlot(Data,sprep=NULL,lys=NULL,plotColors=NULL,
            WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
    data ...
sprep
    ... default is NULL
lys
    ... default is NULL
plotColors
    ... default is NULL
WIDTH
    ... default is 1000
HEIGHT
    ... default is 1414
POINTSIZE
    ... default is 24
MAXARRAY
    ... default is 41

Return
^^^^^^^
.. image:: RawDataSamplePrepControl.png


ratioPlot
++++++++++++++++++++
beta-actin & GAPDH 3'/5' ratio plot

Usage
^^^^^^
ratioPlot(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^^
Data
  ...
quality
  ... default is NULL
experimentFactor
  ... default is NULL
plotColors
  ... default is NULL
legendColors
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^^
.. image:: RawData53ratioPlot_beta-actin.png
.. image:: RawData53ratioPlot_GAPDH.png


RNAdegPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
RNAdegPlot(Data, Data.rnadeg=NULL, plotColors=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
Data.rnadeg
  ... default is NULL
plotColors
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: RawDataRNAdegradation.png

hybridPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
hybridPlot(Data,quality=NULL,plotColors=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
quality
  ... default is NULL
plotColors
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: RawDataSpikeinHybridControl.png

percPresPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
percPresPlot(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
quality
  ... default is NULL
experimentFactor
  ... default is NULL
plotColors
  ... default is NULL
legendColors
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: RawDataPercentPresent.png


PNdistrPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNdistrPlot(Data, WIDTH=1000, HEIGHT=1414, POINTSIZE=24)

Arguments
^^^^^^^^^^
Data
  ...
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24

Return
^^^^^^^
.. image:: RawDataPosNegDistribution.png


backgroundPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
backgroundPlot(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
quality
  ... default is NULL
experimentFactor
  ... default is NULL
plotColors
  ... default is NULL
legendColors
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: RawDataBackground.png


scaleFactPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
scaleFactPlot(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
quality
  ... default is NULL
experimentFactor
  ... default is NULL
plotColors
  ... default is NULL
legendColors
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: RawDataScaleFactors.png


controlPlots
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
controlPlots(Data,plotColors=NULL,experimentFactor=NULL,legendColors=NULL,
  affxplots=TRUE,boxplots=TRUE, WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
plotColors
  ... default is NULL
experimentFactor
  ... default is NULL
legendColors
  ... default is NULL
affxplots
  ... default is TRUE
boxplots
  ... default is TRUE
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: RawDataAFFXControlsProfiles.png


boxplotFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
plotArrayLayout(Data, experimentFactor=NULL, plotColors=NULL, legendColors=NULL, normMeth="",
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41)

Arguments
^^^^^^^^^^
Data
  ...
experimentFactor
  ... default is NULL
plotColors
  ... default is NULL
legendColors
  ... default is NULL
normMeth
  ... default is ""
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24
MAXARRAY
  ... default is 41

Return
^^^^^^^
.. image:: DataBoxplot.png


plotArrayLayout
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
plotArrayLayout(Data,aType=NULL,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24)

Arguments
^^^^^^^^^^
Data
  ...
aType
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24

Return
^^^^^^^
.. image:: RawDataReferenceArrayLayout.png


spatialImages
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
spatialImages(Data, Data.pset=NULL, Resid=TRUE, ResSign=TRUE, Raw=TRUE, Weight=TRUE,
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24)

Arguments
^^^^^^^^^^
Data
  ...
Data.pset
  ... default is NULL
Resid
  ... default is TRUE
ResSign
  ... default is TRUE
Raw
  ... default is TRUE
Weight
  ... default is TRUE
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24

Return
^^^^^^^
.. image:: RawData2DVirtualImage_.png


array.image
++++++++++++++++++++++++++++++++++
alternative function to plot virtual arrays when PLM cannot be run (when > 6 arrays)

Usage
^^^^^^
array.image(Data, pcut=NULL, relative=TRUE, symm=relative,
   balance=relative,quantitative=relative,col.mod=1,postfix="",arrays=NULL,
   WIDTH=1000, HEIGHT=1414, POINTSIZE=24)

Arguments
^^^^^^^^^^
Data
  ...
pcut
  ... default is NULL
relative
  ... default is TRUE
symm
  ... default is relative
balance
  ... default is relative
quantitative
  ... default is relative
col.mod
  ... default is 1
postfix
  ... default is ""
arrays
  ... default is NULL
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24

Return
^^^^^^^
.. image:: RawDataArray.image ... .png


PNposPlot
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot(Data, WIDTH=1000, HEIGHT=1414, POINTSIZE=24)

Arguments
^^^^^^^^^^
Data
  ...
WIDTH
  ... default is 1000
HEIGHT
  ... default is 1414
POINTSIZE
  ... default is 24

Return
^^^^^^^
.. image::

densityFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

densityFunUnsmoothed
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

maFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

nuseFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

rleFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

correlFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

clusterFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::

pcaFun
++++++++++++++++++++++++++++++++++
...

Usage
^^^^^^
PNposPlot()

Arguments
^^^^^^^^^^

Return
^^^^^^^
.. image::


+++++++++++++++++++++++++++++++++++++
Computation & Manipulation Functions
+++++++++++++++++++++++++++++++++++++

getArrayType
++++++++++++++++++++++++++++++++++


addStandardCDFenv
++++++++++++++++++++++++++++++++++

addUpdatedCDFenv
++++++++++++++++++++++++++++++++++

colorsByFactor
++++++++++++++++++++++++++++++++++

deduceSpecies
++++++++++++++++++++++++++++++++++

normalizeData
++++++++++++++++++++++++++++++++++

computePMAtable
++++++++++++++++++++++++++++++++++

createNormDataTable
++++++++++++++++++++++++++++++++++

