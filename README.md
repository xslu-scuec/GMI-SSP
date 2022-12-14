# GMI-SSP

This is a registration method by combining multifeature mutual information with statistical shape prior. It can be faster using ITK multi-threads.

The source code includes two folders: Metrics and Registrations. The users can embed them into elastix platform by using cmake. To be convenient, an executable file based on Windows system can be found from folder 'Release'. 

A bSSFP/T2 cardiac MR registration example is illustrated in folder 'examples'. The "mr_train_060_image.nii.gz" is the fixed image. The "mr_train_060_label.nii.gz" is the fixed mask. The "norm.nii.gz" is the moving image. The other images are the MIND features of fixed and moving images. The "outputpoints.vtk", "meanvector.txt", and "covariance.txt" are statistical shape model files of the left ventricle of the fixed image based on point distribution. 

The users could run:

```
./Release/elastix.exe -f0 ./examples/mr_train_060_image.nii.gz -f1 ./examples/mr_train_060_1.nii.gz -f2 ./examples/mr_train_060_2.nii.gz -f3 ./examples/mr_train_060_4.nii.gz -f4 ./examples/mr_train_060_4.nii.gz -f5 ./examples/mr_train_060_5.nii.gz -f6 ./examples/mr_train_060_6.nii.gz -m0 ./examples/norm.nii.gz -m1 ./examples/norm_1.nii.gz -m2 ./examples/norm_2.nii.gz -m3 ./examples/norm_3.nii.gz -m4 ./examples/norm_4.nii.gz -m5 ./examples/norm_5.nii.gz -m6 ./examples/norm_6.nii.gz -fMask ./examples/mr_train_060_label.nii.gz -fp ./examples/outputpoints.vtk -meanf ./examples/meanvector.txt -covariancef ./examples/covariance.txt -p ./examples/parameters.txt -out ./examples/results
```
