Toolkit for Synaptome Visualization
===================================

Smith Lab @ Stanford University  
<http://smithlab.stanford.edu>

Initial Author: Dan O'Shea [ dan {at} djoshea.com ]  
Source URL: <http://github.com/djoshea/synaptome>


Overview
--------

This code is a collection of related utilities for operating on array 
tomography images of synaptic protein expression in synaptogram format.
A synaptogram is a 3d image stack with multiple named immunohistochemistry
fluorescence channels. Typically these images are displayed in grids where
each channel is a row, and each column represents a serial section, i.e.
one plane lying orthogonal to the z-direction. This code provides means to:

* Load in extracted synaptogram data cubes
* Filter and perform morphological computations on this data
* Calculate features and utilize machine learning to classify synapses
* Visualize feature vectors among different classes in an interactive scatter
  plot
* Display imaging channels and filtered, transformed, transposed, or annotated
  channels in the synaptogram format

The general purpose of this software is to ease the process of working with
data in synaptogram format. It facilitates filtering or manipulating this data
so as to remove noise and locate or quantify morphological structure and using
the results of this filtering and morphological image processing to compute
useful features. These features are easily plugged into supervised learning
algorithms that predict synapse "class" or "type" from their features, which
are ultimately derived from their protein expression characteristics. These
synapses may then be plotted on a scatter plot according to 2 feature
dimensions (although extending this to 3 should be trivial). Each synapse is
a dot colored according to its human provided training label, and optionally
marked with an X if it was misclassified by the last used supervised learning
classifier. Clicking on a synapse dot brings up the synaptogram representation
of that synapse, which is a completely customizable, color annotated version
of the original imaging channels and of filtered channels computed from the
imaging data, from which the feature vectors are ultimately extracted. This
workflow allows for new feature extraction schemes to be implemented and
evaluated quickly and for error cases and outliers to be studied in detail.
This process also allows for dataset exploration, in which filtered and
annotated channels highlight the structure of the data and provide more
helpful information than raw synaptograms.

Setting Up the Path and the 'ds' Structure
------------------------------------------

Most of the code lies in the root synaptome directory, whereas code specific
to a particular dataset, such as a particular loading script or filtering
script, is generally located in a directory used only for that dataset. For
this reason, you should *add the synaptome directory with subdirectories (or
at least its filters subdirectory) to your matlab path* and then operate out 
of the local dataset directory. 

Synaptome makes heavy use of Matlab structures to store and organize data. In
particular, a variety of information related to a given dataset is stored
within a structure referred to as <tt>ds</tt> within the code. Most functions take 
<tt>ds</tt> as their first argument, and will often write data into the <tt>ds</tt>
structure and pass along the augmented ds structure as the return value, allowing 
these functions to be called as:

    ds = somefunction(ds, otherargs, ...);

Some key fields of the <tt>ds</tt> structure, organized by category, are described
below.

Example Usage
-------------

Add the synaptome directory and all subdirectories to the Matlab path.

Load in a preloaded dataset as a .mat file.

	cd /path/to/synaptome/kdm200
	load trainds.mat

Run the example filtering and feature extraction script.

	features_kdm200watershed

Click on a few points in the scatter plots to view the annotated synaptograms.

Run a simple supervised learning algorithm to evaluate the utility of the
features computed for classification.

	testft

Imaging Data
------------

Within imaging data arrays, the dimensions are ordered as [Z Y X], which
makes, image view orientation of pixels identical to the standard matrix
orientation (i.e. Z addresses pages, Y addresses rows, X address columnds)

### Fields

<b><tt>dsname</tt></b> : name of dataset or imaging session

<b><tt>chlist</tt></b> : names of original imaging channels present in the dataset

<b><tt>nimch</tt></b> : number of original imaging channels

<b><tt>sgdim</tt></b> : [Z Y X] dimensions of synaptogram

<b><tt>sgaspect</tt></b>: [Z Y X] aspect ratio of each image voxel

<b><tt>sg</tt></b> : struct with individual synapse structures beneath

<b><tt>sg.im</tt></b> : image channel data (4d: nimch x sgdim)

<b><tt>sg.str</tt></b> : descriptive string for this synapse (coordinatess, label, etc.)

<b><tt>sg.im</tt></b> : image channel data (4d: nimch x sgdim)

<b><tt>sg.str</tt></b> : descriptive string for this synapse (coordinatess, label, etc.)

### Functions

#### <tt>chdat = getimchannel(ds, chname)</tt>

Returns data from a specific imaging channel by name for all synapses
chname is the name of the original imaging channel to search in chlist.

Filtered and Color Annotated Channels
-------------------------------------

Filtered imaging data resides in ds.trd (as in "training data") instead of
<tt>ds.sg(:).im</tt>. This is grayscale data of the same size (<tt>ds.sgdim</tt>) and aspect ratio 
(<tt>ds.sgaspect</tt>) as the original imaging channels. Typically, the original imaging
data are copied verbatim to the filtered channel array (<tt>ds.trd</tt>). Filtered
channels are used for both feature computation and for visualization.

Color annotation channels are essentially red green blue versions of filtered
channels. They are used solely for visualization as rows in the synaptograms.
For example, once distinct puncta of a certain protein are isolated in a
specific channel, they can be colored so as to distinguish each puncta from
all others according to a colormap. The result of this coloring would be saved
as a color annotation channel and referred to by name when building the
synaptogram.

Both types of channels are typically generated by using the filtdat function
and added to ds using the addchannel or addcolorchannel function, which are
described in detail below. There are a variety of filters within filtdat.m. To
get a sense of how to effectively utilize and chain filters in filtdat.m, see
one of the sample <tt>features_kdm200.m</tt> or <tt>features_gw300.m</tt> scripts 
along with the detailed documentation within filtdat.m adjacent to each filter.

### Fields

<b><tt>trd</tt></b> : number of data channels x ntrain x sgdim (5d array)

<b><tt>trdlist</tt></b> : names of each data channel (same as first dimension of trd)

<b><tt>colorch</tt></b> : number of color channels x ntrain x sgdim x 3 (6d array last is RGB)

<b><tt>colorchname</tt></b> : names of color channels (same as first dimension of colorch)

### Functions

#### <tt>ds = addchannel(ds, dat, dname, normalize)</tt>

Adds a filtered image channel to ds.trd with name dname.

<tt>dat</tt> : ntrain x sgdim data to store as channel

<tt>dname</tt> : name string for new channel

<tt>normalize</tt> : boolean, indicates whether all data should be as a whole linearly shifted to fit the [0, 1] range. defaults to 1.

#### <tt>ds = addcolorchannel(ds, dat, dname)</tt>

Adds a color channel to ds.colorch with name dname. 

<tt>dat</tt> : ntrain x sgdim x 3 RGB color values, all values in range [0, 1].

#### <tt>[filt info] = filtdat(ds, dat, params)</tt>
	
Performs some kind of filtering or computational operation on some
imaging data and returns the results and optionally some metadata.
filtdat merely returns the filtered output and metadata, it is not 
designed to make any changes directly to ds. This enables filters to be
chained easily before adding the output to the trd list via addchannel.
The filtered output may be easily added back to ds.trd via:

        filtparams = struct('type', 'typename', 'option1', value1);
        ds = addchannel(ds, filtdat(ds, 'InputChannelName', filtparams), ...
                'OutputChannelName');

Individual filters are typically implemented as typename_filt.m files in
the filters directory (but may be placed anywhere on the matlab path).
They may also be implemented within the filtdat.m where indicated.
	
See the comments for each individual filter function within filtdat.m for
more specific documentation on usage.

<tt>dat</tt> : the primary input to the operation, specified either directly
as an array ntrain x sgdim or by name (to be retrieved using getchannel

<tt>params</tt> : a struct specifying the details of the operation to perform
params.type specifies which operation to perform and must be 
specified. Other values under params depend on the operation being
performed. These values can include parameters used for all synapses
globally, arrays of parameters used for individual synapses, or even
entire additional imaging data arrays.

<tt>filt</tt> : the primary output of the operation as an array ntrain x sgdim

<tt>info</tt> : a structure containing metadata returned by some operations.
These fields may contain arrays of values generated for each synapse
individually or secondary output channel arrays.

#### <tt>ds = getchannel(ds, dname, syn)</tt>

Retrieves data for a particular channel for all or a particular
synapse. Searches trdlist first (filtered channels), then original
imaging channels (using getimchannel), and then colorchannels. Error
if channel still not found.

<tt>dname</tt> : name of the channel for which to search.

<tt>syn</tt> : optional id number in the interval [1 ds.ntrain]. If used,
the returned data will be for the specified synapse only, instead of
an array for all synapses.

Visualization
-------------

Individual synapses may be visualized as a synaptogram, a 2D grid format where
each row is a channel and each column is a serial section. The rows may
represent either filtered channels (typically image channels are copied over
into filtered channels, essentially transferring them from ds.sg.im to ds.trd
for convenience) or color annotated channels. 

The rows used to build a synaptogram for any of the synapses are specified in
the ds.vis cell array, where each element is a channel name or an array of up
to 3 channel names which will be merged as RGB to make a color composite row on
the fly. Channels which are already color annotated may only be used by
themselves and not included in an RGB composite row. Typically visualizations 
are constructed by clearing the ds.vis and ds.visname lists and using addvisrow
to add each row. For example, to view Synapsin in row 1 as grayscale, VGlut1 as 
red and VGlut2 as blue in row 2, and a color annotation channel Synapsin_puncta 
in row 3, use:

    ds.vis = {};
    ds.visname = {};
    ds = addvisrow(ds, 'Synapsin');
    ds = addvisrow(ds, {'VGlut1', 'VGlut2'});
    ds = addvisrow(ds, 'Synapsin_puncta');

### Fields

<b><tt>vis</tt></b> : array of channel names (or arrays of 3 channel names for RGB composites)
used for each row of the synaptogram visualization

<b><tt>visname</tt></b> : array of names printed next to each row of the synaptogram.
Typically these are automatically filled out when addvisrow() is used.

### Functions

#### <tt>ds = addvisrow(ds, visdat, name)</tt>

Used to add a row to the ds.vis list.
	
<tt>visdat</tt> : either the name of a channel (filtered or color annotation)
or an array of 3 channel names. If one of the colors is not
specified (i.e. array of length 2) or is given by '', the color will
be set to 0 in the compositing.

<tt>name</tt> : an optional string used to override the default name printed
next to each row on the synaptogram. The default is the original
name or each of the channel names joined by ' / ' in between.

#### <tt>syntable(ds, i)</tt>

Displays a data table containing the value of each feature for a
particular synapse, typically evoked within viewsg to be displayed
alongside the synaptogram.

<tt>i</tt> : synapse index in the interval [1 ds.ntrain]

#### <tt>viewsg(ds, i, fignum)</tt>

Displays the synaptogram for synapse i in figure(fignum). This is the
method called from a synplot's on click callback function to bringup
the synaptogram for the clicked point. Uses ds.vis and ds.visname to
build up the synaptogram image row by row using viewsgrow and then to
display the full synaptogram on a dark background with text labels to
either side.

<tt>i</tt> : synapse index in the interval [1 ds.ntrain]

<tt>fignum</tt> : number of figure to display in. default is gcf.

#### <tt>rowdat = viewsgrow(dat, bordercol, vspacing, hspacing)</tt>
	
<tt>dat</tt> : either a 3d intensity array in [z y x] coords, or a
struct where dat.R, dat.G, dat.B are 3d intensity arrays in [z y x]
coords (for RGB color)

<tt>bordercol</tt> : 3x1 RGB value used to fill the margins around each XY tile

<tt>vspacing</tt> : padding above and below each tile

<tt>hspacing</tt> : padding left of each tile and to the right of the last tile

Supervised Learning
-------------------

The assignment of synapses into classes is complicated because of variability
in whether synapses are allowed to belong to multiple classes. The distinction
is delineated in the code as "classes" versus "labels". A neuron can be
assigned to multiple classes, e.g. glut1 and glut2, but it belongs to one and
only one label, e.g. "glutBoth". Within the kdm200/loadds.m sample loading
script, there are 3 classes and 5 labels, and a special labelmatch field
contains a 3-dim vector for each of the labels indicating which combination of
classes each label represents. This approach is confusing, but resulted from
the different ways in which human voting was performed. The classes list
represents human names for different known types of synapses, and some
synapses ambiguously seemed to express markers indicative of multiple classes.
To make the supervised learning problem straightforward, labels were
introduced such that each neuron belonged to exactly 1 label. However, this
system must be taken into account when switching from an n-class machine
learning problem to a 2 class problem such as glut1 vs. NOT glut1, which must
properly assign glut1 membership to neurons marked as glutBoth.

See testft.m for an example of training a simple classifier on the training
set and evaluating it's training and cross-validation error.

### Fields

<b><tt>classes</tt></b> : names for each class

<b><tt>labelnames</tt></b> : names for each label

<b><tt>labelcolors</tt></b> : RGB values used to represent each label

<b><tt>ntrain</tt></b> : number of synapses in training set

<b><tt>votes</tt></b> : full voting table by classes (synapses x classes x voters)

<b><tt>nvotes</tt></b> : number of voters

<b><tt>votesbylabel</tt></b> : voting table by labels (synapses x labels x voters)

<b><tt>trainlabel</tt></b> : actual label id assigned to each synapse (indexes labelnames)

<b><tt>trainlabelrunnerup</tt></b> : label id with the second most votes

<b><tt>trainlabelconf</tt></b> : confidence measure determined by (votes to
winner - votes to runner up) / nvotes

<b><tt>trainactive</tt></b> : boolean indicator of whether to include this synapse in
supervised learning training and testing, typically this is a thresholded
version of trainlabelconf to take only examples upon which humans agree

### Functions

#### <tt>testft</tt>

Sample script for performing cross-validation error tests using a simple linear classifier.


Feature Computation
-------------------

For each synapse in the training set, a feature vector may be computed for use
in supervised learning. Each feature is a scalar value typically computed
based on that channels filtered image data. One example is integrated
brightness, which is simply the sum of all voxel values in a particular
channel for that synapse. Each feature has a name, which is used throughout
the code to refer to and retrieve the value of that feature for particular
synapses.

### Fields

<b><tt>ft</tt></b> : ntrain x number of features, rowwise matrix of feature vectors

<b><tt>ftname</tt></b> : list of feature names

### Functions

#### <tt>ds = addfeature(ds, ftdat, name)</tt>

Adds a new feature columnn to ds

<tt>ftdat</tt> : array of scalar values for the feature column for each
of the ntrain synapses

<tt>name</tt> : unique string used to refer to that feature

#### <tt>[ftmat idx] = getfeature(ds, names)</tt>

Retrieves a vector or matrix of values for each synapse for one or
more features by name.

<tt>names</tt> : either a string or an array of strings

<tt>ftmat</tt> : either a vector or matrix of feature values for all synapses

<tt>idx</tt> : column index into the feature matrix <tt>ds.ft</tt> at which each
named feature was found

Feature Scatter Plots
---------------------

Once features vectors have been computed for all training synapses, the
training set may be visualized as a 2D scatter plot where each synapse is
represented as a dot colored according to its human provided trainlabel. The
two axes are specified by named features within ds.ft. 3D plots are
not implemented but should be trivial if the need arises. Optionally, synapses
which were misclassified by the most recently trained supervised learning
classifier may be marked with an X. Clicking on a synapse dot in this plot
will load a synaptogram visualization for that particular synapse.

### Functions

#### <tt>synplot(ds, features, errors, twoclass)</tt>

Displays the interactive scatter plot of feature values
features is an array of strings representing feature names to be used
for the X and Y axes. 
	
<tt>errors</tt> : binary vector indicating which synapses to mark with an X, 
defaulting to <tt>zeros(ds.ntrain,1)</tt>

<tt>twoclass</tt> : array of label names used when the scatter plot should 
color synapses according to only 2 categories, those whose label is 
in twoclass (as red) and all others (as black)

