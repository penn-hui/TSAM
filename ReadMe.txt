A simple introduction:

TSAM is the abbreviation of Two-Step-Averaging-Method which implements prediction or classification
with a two step process. In the first step, a xgboost model is adopted to predict the first-step
scores. Then, the xgboost selected important features will be combined with our profiled Hidden
Markov Model features to predict the second-step score. Finally, the TSAM output a averaged score
as the predicted efficiency score.

*************************************************************************************************************
*************************************************************************************************************
What's included in the tool package:

The offline TSAM includes the Python version and the Matlab version. Both of these two version
tool can be used to predict sgRNAs' cutting efficiencies or classify the sgRNAs into high or
efficent ones. We did not control the precision of the floating-point operations, thus the predicted
scores of the Matlab version and the Python version may contain some differences (within +/-0.01).
The Matlab version can just be used to do the prediction and classfication but can not used to implement
the cross-validation and the independent test that introduces in our journal paper. The Python version
contains the function to implement those experiments. More details can be found in the ReadMe.txt file
in the folder python_codes.

*************************************************************************************************************
*************************************************************************************************************
Package dependency:

The following packages must be installed to run the Python and Matlab version tool:
1. xgboost
2. biopython
3. numpy
4. sklearn
5. libsvm
6. scipy

More details about how to install these packages can be found in the ReadMe.txt file in the matlab_codes and
python_codes folder.

*************************************************************************************************************
*************************************************************************************************************
How to prepare the inputs:
 
We have described the steps to prepare the input sequences in the ReadMe.docx file.

*************************************************************************************************************
*************************************************************************************************************
How to run the codes:

1. How to run the Python version codes:

After downloading the codes and decompressing it, one can open a console application (command prompt on the
windows os or a console on the linux os). Using the python_codes folder path as the current working path and
run the codes. For example:
     
     cd python_codes
     python TSAM_python.py ../example_input_files/test1.fa annotated 1 1 1 1   

(please pay attention on the path character ‘/’ on linux os and ‘\’ on windows os)

Then, you can find the .csv file in the folder: python_codes/predicted_scores/ predict_results.csv
More details about how to run the python codes can be found in the ReadMe.txt file under the python_codes folder 
********************************************************************************************************************
2. How to run the Matlab version codes:

Download the codes and decompress the files, then open a Matlab software (Matlab 2015 or higher). Adding the folder
TSAM/matlab_codes/ into the working path. Then run the prediction command. For example:
     
      Predict_score=TSAM(../example_input_files/test1.fa, 1, 1, 1);   

(please pay attention on the path character ‘/’ on linux os and ‘\’ on windows os)
More details about how to run the python codes can be found in the ReadMe.txt file under the matlab_codes folder 

********************************************************************************************************************
********************************************************************************************************************

***Please contact Hui Peng: Hui.Peng-2@student.uts.edu.au if you encounter some problems when run the codes.