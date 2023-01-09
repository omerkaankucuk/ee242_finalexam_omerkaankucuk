outputFolder=fullfile('caltech101');
rootFolder= fullfile(outputFolder,'101_ObjectCategories');%103 categories

categories = {'leaf_seg','airplanes', 'ferry' ,'laptop','theo_out','real_out'};

imds = imageDatastore(fullfile(rootFolder,categories),'LabelSource','foldernames');
tbl = countEachLabel(imds);%label & numberof data
minsetcount = min(tbl{:,2});%choose min number of data 

imds=splitEachLabel(imds,minsetcount,'randomize');%randomly choose 67 sample each class
countEachLabel(imds)

leaf = find(imds.Labels == 'leaf_seg',1);
airplanes = find(imds.Labels == 'airplanes',1);
ferry = find(imds.Labels == 'ferry',1);
laptop = find(imds.Labels == 'laptop',1);

theoout=find(imds.Labels == 'theo_out',1);
realout=find(imds.Labels == 'real_out',1);
% custom training with resnet in matlab
figure
subplot(2,2,1);
imshow(readimage(imds,leaf));
subplot(2,2,2);
imshow(readimage(imds,airplanes));
subplot(2,2,3);
imshow(readimage(imds,theoout));
subplot(2,2,4);
imshow(readimage(imds,realout));



net=resnet50();
figure
plot(net)
title('Architecture of ResNet')
set(gca,'Ylim',[150 170]);

net.Layers(1)% inputsize [224 224 3]
net.Layers(end)%trained 1000 times 

numel(net.Layers(end).ClassNames)
[trainingSet,testSet]= splitEachLabel(imds,0.3,'randomized');
%for traning %30 validation %70

%Before sending to Cnn we need to modify inputs
imageSize=net.Layers(1).InputSize;
%first perform on the training data set
augmentedTrainingSet=augmentedImageDatastore(imageSize,...
    trainingSet,'ColorPreprocessing','gray2rgb');
%second perform on the test data set
augmentedTestSet=augmentedImageDatastore(imageSize,...
    testSet,'ColorPreprocessing','gray2rgb');



w1=net.Layers(2).Weights;
w1= mat2gray(w1);%matrix to grayscale image
figure
montage(w1)
title('First Convolutional Layer Weight')
%these features are processed by deep network layers which combine
%the early features to form the higher level image features
%these high-level features are better for recognition tasks
%because they combine all the primitive features into a
%better image representation using activition method we can easly
%extract features one of the deeper layers we can choose any layer
%however I prefer extracting features from the layer right before the
%classification layer in ResNet-50 this layer named (FC 1000)

%lets extract hte features
featureLayer = 'fc1000';

%use activation method(ensure that the CNN and the image data fit into GPU
%first GPU if not CPU 
%output arranged columnar output speeds up SWM training
trainingFeatures = activations(net,...
    augmentedTrainingSet,featureLayer,'MiniBatchSize',32,'OutputAs','columns');
%now we need the levels of traing set
trainingLables = trainingSet.Labels;

%Fit-Class Errror Correcting Output Codes (ECOC)
% Returns full trained multi class error correcting output
% K(K-1)/2 binary support vector machine
% K = number of unique class names
%returns a trained model
classifier=fitcecoc(trainingFeatures,trainingLables,...
    'Learner','Linear','Coding','onevsall','ObservationsIn','columns');
%extract the features from test
%using these test feature we can measure the accuracy of the
%trained classifier if we pass the test features to classifier and compare
% and compare it with obtained features we can calculate the accuracy of
% classifier to do it we use predict func
testFeatures = activations(net,...
    augmentedTestSet,featureLayer,'MiniBatchSize',32,'OutputAs','columns');

%predict func returns a vector of predicted class levels based on the
%trained classifier

%predicted labels
predictLabels = predict(classifier,testFeatures,'ObservationsIn','columns');
%real labels
testLables = testSet.Labels;

%confusion matrix to evaluate the performance of the classifier
%first argument is the known value(testLabels,predictLabels)
confMat=confusionmat(testLables,predictLabels);

%it converted value of confMat in to percantage(1 = %100,0 = %0)
confMat1=bsxfun(@rdivide,confMat,sum(confMat,2));


%calculated the entire row and put the vaalue as first element 
%first element is 47 ,47,47 and do it for every col
%sum(confMat,2)

mean(diag(confMat1)) % ans should be 0.xxxx




%*********** Classify any image **************

newImage = imread(fullfile('test.png'));

ds = augmentedImageDatastore(imageSize,...
    newImage,'ColorPreprocessing','gray2rgb');

newimagefutures = activations(net,...
    ds,featureLayer,'MiniBatchSize',32,'OutputAs','columns');

label1 = predict(classifier,newimagefutures,'ObservationsIn','columns');

sprintf('The loaded image belongs to %s class',label1)

%******* TEST AREA *********

% [C,scores] = semanticseg(newImage,net);
% 
% B = labeloverlay(newImage, C);
% figure
% imshow(B)




% semanticseg
%https://www.mathworks.com/help/vision/ref/semanticseg.html?s_tid=doc_ta#d124e183473
%https://www.mathworks.com/help/releases/R2020b/deeplearning/ug/create-simple-semantic-segmentation-network-in-deep-network-designer.html

%https://www.mathworks.com/help/deeplearning/ref/resnet50.html
%https://www.mathworks.com/help/deeplearning/ref/trainnetwork.html


































