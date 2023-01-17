aprox_hermit: main.o splines.o points.o hermit.o rozwiazywacz.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o hermit.o rozwiazywacz.o -L gaus -l ge

aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c
hermit.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c hermit.c

rozwiazywacz.o: rozwiazywacz.c rozwiazywacz.h
test1: aprox_hermit clean
	./aprox -s spl -p test/dane.1 -g myplot -f 5.1 -t 5.7 -n 300
test2: aprox_hermit clean
	./aprox -s spl -p test/dane.2 -g myplot -f 0 -t 1 -n 300
test3: aprox_hermit clean
	./aprox -s spl -p test/dane.3 -g myplot -f 0 -t 1 -n 300
	
test4: aprox_hermit clean
	./aprox -s spl -p test/dane.4 -g myplot -f 1 -t 2 -n 300
	
test5: aprox_hermit clean
	./aprox -s spl -p test/dane.5 -g myplot -f 0 -t 1 -n 500
	
test6: aprox_hermit clean
	./aprox -s spl -p test/dane.6 -g myplot -f 5 -t 9 -n 300
	

	
	
	

.PHONY: clean

clean:
	-rm *.o 
