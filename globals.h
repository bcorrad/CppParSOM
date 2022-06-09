#ifndef UNTITLED_GLOBALS_H
#define UNTITLED_GLOBALS_H

constexpr int N_THREADS = 1; // [1, 2, 3, 4].
constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;
constexpr int WEIGHT_SIZE = 100; // size of weight vector of each node [10, 100, 1000].
constexpr int EPOCHS = 100;   //the number of epochs desired for the training [10, 100, 1000].
constexpr double START_LR = 0.001;    //the value of the learning rate at the start of training.
constexpr int GRID_ROWS = 10;
constexpr int GRID_COLS = 10;
constexpr int N_INPUTS = 1000; // number of weights each node must contain. [100, 1000, 10000]

#endif //UNTITLED_GLOBALS_H
