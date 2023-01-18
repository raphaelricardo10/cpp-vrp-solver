lib:
	g++ --shared -o ga.so routes.cpp -fPIC

lib-mac:
	clang++ -std=c++20 -stdlib=libc++ --shared routes.cpp -o ga.so

main-mac:
	clang++ -std=c++20 -stdlib=libc++ main.cpp -o main