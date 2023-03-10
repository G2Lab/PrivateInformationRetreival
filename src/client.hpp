#pragma once

#include <iostream>
#include <helib/helib.h>

class Client{
public:
    Client(const helib::Context &context);
private:
    const helib::Context* context;
};