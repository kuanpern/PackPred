CC = g++
CCFLAGS = -lstdc++ -Wall -Wextra  
bindir = /usr/local/bin/
libdir = /usr/local/lib/depth-2.0/
lib = par/ water_models/
executables = DEPTH complex-depth compute_deviation predict-binding-site generate_sites ASA g_distmatrix
INSTALL_H = src/install_dir.h

all:
	echo "install program in current directory"
	echo "#define EXE_DIR string(\""$(PWD)"\")" > $(INSTALL_H)
	echo "#define LIB_DIR string(\""$(PWD)"\")" >> $(INSTALL_H)
	$(CC) $(CCFLAGS) -fopenmp src/main.cpp -o DEPTH
	$(CC) $(CCFLAGS) src/complex-depth.cpp -o complex-depth
	$(CC) $(CCFLAGS) src/predict-binding-site.cpp -o predict-binding-site
	$(CC) $(CCFLAGS) src/ASA.cpp -o ASA
	$(CC) $(CCFLAGS) src/g_distmatrix.cpp -o g_distmatrix

install:
	echo "install program in " $(bindir)
	echo "#define EXE_DIR string(\""$(bindir)"\")" > $(INSTALL_H)
	echo "#define LIB_DIR string(\""$(libdir)"\")" >> $(INSTALL_H)
	$(CC) $(CCFLAGS) -fopenmp src/main.cpp -o DEPTH
	$(CC) $(CCFLAGS) src/complex-depth.cpp -o complex-depth
	$(CC) $(CCFLAGS) src/predict-binding-site.cpp -o predict-binding-site
	$(CC) $(CCFLAGS) src/ASA.cpp -o ASA
	$(CC) $(CCFLAGS) src/g_distmatrix.cpp -o g_distmatrix
	rm -rf $(libdir); mkdir $(libdir); cp -r $(lib) $(libdir)
	cp $(executables) $(bindir)

clean:
	rm $(executables); rm $(INSTALL_H)

remove:
	cd $(bindir); rm $(executables) 
	rm -rf $(libdir)

