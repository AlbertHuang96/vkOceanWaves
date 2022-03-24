# vkOceanWaves

Vulkan implementation of ocean waves in both DFT and FFT

It is now basically a well organized version of Discrete Fourier Transformation Ocean Wave. 
The base project (base.lib) from the SaschaWillems examples codebase was used as a framework to provide a wrapper of Vulkan APIs.

It is very slowly due to the algorithm - discrete Fourier Transform.

In the latest version I update it to the Fast Fourier Transform algorithm. 
Under NVIDIA Geforce RTX 3050 Ti Laptop GPU and Intel i7-11800H @2.30GHz, the FPS is about 120+ in debug mode.

Updated: In release mode it is about 300 fps.

You can use mouse and wheel to navigate. For now I do not set any cull mode.

NOTICE: 
You can select your GPU in the line 870 in vulkanexamplebase.cpp. 
The variable selectedDevice is 1 for default. 
You can change it to 0 in case:
 you want to switch to another GPU 
 or you only have only one GPU and the program crashes at the GPU selection.
 (The GPU will show on the window title and the control panel as well.)
 
![wave](https://user-images.githubusercontent.com/42757470/159940977-d4095f51-8ebb-47dd-9f47-5e612a761801.png)
