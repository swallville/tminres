serial :
	cd SerialExample; $(MAKE) all;
	
trilinos :
	cd TrilinosExample; $(MAKE) all;

dox :
	doxygen Doxyfile.in;
	
clean :
	cd SerialExample; $(MAKE) clean;
	cd TrilinosExample; $(MAKE) clean;
	
clean-dox :
	rm -r doc;


RELEASE_DIR = tminres-0.1
release :
	$(MAKE) dox;
	mkdir $(RELEASE_DIR);
	cp *.hpp $(RELEASE_DIR);
	cp README $(RELEASE_DIR);
	cp Doxyfile.in  $(RELEASE_DIR);
	cp Mainpage.dox $(RELEASE_DIR);
	cp Makefile $(RELEASE_DIR);
	mkdir $(RELEASE_DIR)/SerialExample;
	cp SerialExample/*.hpp $(RELEASE_DIR)/SerialExample;
	cp SerialExample/*.cpp $(RELEASE_DIR)/SerialExample;
	cp SerialExample/Makefile $(RELEASE_DIR)/SerialExample;
	mkdir $(RELEASE_DIR)/TrilinosExample;
	cp TrilinosExample/*.hpp $(RELEASE_DIR)/TrilinosExample;
	cp TrilinosExample/*.cpp $(RELEASE_DIR)/TrilinosExample;
	cp TrilinosExample/*.mtx $(RELEASE_DIR)/TrilinosExample;
	cp TrilinosExample/Makefile $(RELEASE_DIR)/TrilinosExample;
	cp TrilinosExample/Makefile.in $(RELEASE_DIR)/TrilinosExample;
	cp -r doc $(RELEASE_DIR)/doc;
	rm -rf $(RELEASE_DIR)/doc/html/formula.repository;
