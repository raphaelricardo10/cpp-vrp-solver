libpath = src/lib/ga

lib:
	g++ --shared -o ${libpath}/ga.so ${libpath}/routes.cpp -fPIC

lib-mac:
	clang++ -std=c++20 -stdlib=libc++ --shared ${libpath}/routes.cpp -o ${libpath}/ga.so

main-mac:
	clang++ -std=c++20 -stdlib=libc++ ${libpath}/main.cpp -o ${libpath}/main