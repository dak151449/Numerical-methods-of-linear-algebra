#include "Func_image.h"

void Func_image::plot_image(std::vector<std::vector<double>> image) {
    if (image.size() != 2) {
        std::cout << "Error: image must be 2D" << std::endl;
        return;
    }
    std::vector<double> x_axis = image[0];
    std::vector<double> y_axis = image[1];
    plt::plot(x_axis, y_axis);
    plt::show();
}

void Func_image::plot_image(func F) {
    
    std::vector<double> x_axis;
    std::vector<double> y_axis;

    for (float i = -100; i < 100; i += 0.1)
    {
        x_axis.push_back(i);
        y_axis.push_back(F(i));
    }
    
    plt::plot(x_axis, y_axis);
    plt::show();
}
