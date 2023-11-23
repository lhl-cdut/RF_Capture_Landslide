All code is run on Matlab.

If You Have Your Own Dataset:

1. Preprocess your data initially.
2. Perform feature extraction using "data_processing (dataset)".
3. Save the final attributes as a ".csv" file.
4. Train the model using "model_training". The resulting model will be saved as a ".mat" file.
5. If you wish to make predictions on new data, open "newdatapredict". This includes signal processing, feature extraction, and loading your own recognition model.

If You Don't Have Your Own Dataset:

We provide simple models in "Example_models" and pre-extracted feature data in "feature_20," available in "model_training."
You can also directly classify new data in "newdatapredict."
Keep in mind that the provided models are for code validation purposes. If you have specific requirements, feel free to adjust the code accordingly.