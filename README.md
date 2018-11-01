## Packpred
Predicting the functional consequences of point mutations in proteins.

### Server demo
 - A running server (currently supported by IISER Pune) can be found at http://cospi.iiserpune.ac.in/packpred/

### References
```
Status: submitted
Title: Packpred: Prediction the functional strengths of point mutations in proteins.
Authors: Kuan Pern Tan, Parichit Sharma, Kwoh Chee Keong, M.S. Madhusudhan
```

### Installations
Assume OS Ubuntu 16.04
```
# install dependencies
$ sudo apt-get install python2.7 g++
```

#### Virtual environment
```
$ cd $install_dir
$ virtualenv -ppython2.7 venv
$ venv/bin/pip install -r requirements.txt
```

#### Modeller
Modeller program needs to be installed on the system. Please refer to https://salilab.org/modeller/download_installation.html for installation guide for different opearating system.

```
# compile the dependency programs
$ cd $install_dir/exec/depth-2.0/bin/; make
$ cd $install_dir/exec/YETI-2.0/; python2.7 install.py

# link with modeller libraries
$ cd $install_dir/venv/lib/python2.7/site-packages
$ ln -sf /usr/lib/python2.7/dist-packages/modeller/ .
$ ln -sf /usr/lib/python2.7/dist-packages/_modeller.so .
```

### test-run
```
$ cd $install_dir/tests/unit_test_b4/
$ $install_dir/venv/bin/python $install_dir/packpred_runner.py --input input.json --exedir $install_dir --workdir $PWD
```
