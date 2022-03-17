Script for the conactenated autoencoder (AE)
---------------------------------------
* The concatenated AE is used for the comparison with scAEGAN. 
* Consists of two encoders for each respective domains, which are then projected down to the bottleneck layer.
* The bottleneck layer contains the integrated low-dimensional representation of the two domains.
* To run the the the concatenated AE script run the following command and provide the dataset paths in the AE_concatenated.py
```
python AE_concatenated.py
```
