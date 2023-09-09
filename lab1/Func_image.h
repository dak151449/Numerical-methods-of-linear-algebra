#pragma once
#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;


typedef double (*func)(double a);


class Func_image {
    public:
        Func_image(){};
        static void plot_image(std::vector<std::vector<double>> image);
        static void plot_image(func F);
};