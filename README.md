# NARU_PARREC2NIFTI
An SPM plugin to handle converting Philips PAR/REC files to NIFTI, with support for Dual Echo protocols.

## Installation
Just take all the Matlab code in this repository and stick it into your spm/toolbox directory. For example, on my system, the code is installed to `C:\Users\mbmhscc4\MATLAB\Toolboxes\spm12\toolbox\PARREC2NIFTI`. Everything prior to `spm12/` will of course be different for you. This has not been tested with SPM8, but I do not forsee compatibility issues.

If you do not use Git, the easiest way to get the code is to download everything as a .zip file. The option to do so is near the top of this page: click "Clone or download" and select "download zip".

![downloadsource](https://user-images.githubusercontent.com/2721307/27448839-2edf8abc-577e-11e7-91cd-ae5baa57f314.PNG)

## How to find it
If you've installed it properly, you should be able to access this tool through the SPM `batch` interface. Open the `batch` interface and select "SPM>Tools>Convert PAR to NIFTI" from the menu.

![howtofindit](https://user-images.githubusercontent.com/2721307/27447458-f7b981e6-5778-11e7-939e-9945fdbf6dac.PNG)

## How to use it
The interface will initially give you two options: to add a dataset and specify an output directory. A dataset is a collection of PAR/REC files with a similar acquisition protocol (EPI, Dual Echo EPI, B0, or T1), and which you want output NIFTI versions into a common directory. Each PAR/REC conversion is independent of all the others, so this grouping of PAR/REC pairs into datasets is really just a way to process several files through the same conversion procedure.

![initialdisplay](https://user-images.githubusercontent.com/2721307/27447460-f7bbec2e-5778-11e7-8f4a-6c05952580e9.PNG)

## Chosing a scan protocol
To see the available options, click on `Scan Protocol`. Select from the window below.

![choosescanprotocol](https://user-images.githubusercontent.com/2721307/27447461-f7bd7210-5778-11e7-87bb-0f7130fe36b2.PNG)

### Dual Echo
If the scans were acquired with a dual echo protocol, make sure to select that from the set of options that appear once you highlight `Scan Protocol`.

![dualechoselection](https://user-images.githubusercontent.com/2721307/27447459-f7bba854-5778-11e7-96c2-adc04f29f0e8.PNG)

Chances are, you have a short and a long echo, and the default labels work well for you. If for some reason you want to change them, or you used some crazy multi-echo protocol, you can change the labels and even add more. The important thing is that the number of labels you provide is the same as the number of echoes in your protocol. If there is a mismatch, the process will abort with an error.

A second important thing to keep in mind when editing the echo labels is that they should be provided in order of ascending length. It is relevant that 'short' precedes 'long' in the default definition. Short echos will always be processed before longer echoes.

![dualecholabelentry](https://user-images.githubusercontent.com/2721307/27447462-f7bed754-5778-11e7-8c15-611bcfe87e2e.PNG)

### B0 (phase/magnitude map)
Files of this kind will two have different scan types, one that provides an estimate of magnitude and a second that provides an estimate of phase. The conversion function will identify these automatically and produce two NIFTI files that are appropriately labeled. These labels (phase and magnitude) cannot be altered through this interface.

## load_parrec.m
The function `load_parrec()` is the only bit of the code that actually interacts with the PAR/REC files. If you have questions or concerns about how information is being extracted from these files, then please consult that code and associated comments. Or, of course, you know where I sit and how to contact me.
