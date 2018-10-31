## Packpred
Predicting the functional consequences of point mutations in proteins.

### References
```
Status: submitted
Title: Packpred: Prediction the functional strengths of point mutations in proteins.
Authors: Kuan Pern Tan, Parichit Sharma, Kwoh Chee Keong, M.S. Madhusudhan
```

### installations
Modeller program needs to be installed on the system. Please refer to https://salilab.org/modeller/download_installation.html for installation guide for different opearating system.

Assume OS Ubuntu 16.04
```
# install dependencies
$ sudo apt-get install python2.7 g++

# compile the dependency programs
$ cd $install_dir/exec/depth-2.0/bin/; make
$ cd $install_dir/exec/YETI-2.0/; python2.7 install.py
```

### test-run
```
$ cd $install_dir/tests/
$ $install_dir/venv/bin/python $install_dir/packpred_runner.py --input input.json
```

### Server instances
 - A running server supported by IISER Pune can be found at http://cospi.iiserpune.ac.in/packpred/
