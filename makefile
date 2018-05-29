all:
	cd src; make

debug:
	cd src; make debug

clean:
	rm -f src/*.o src/*.mod bin/LES2
