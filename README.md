# vkOceanWaves
---
#### Description  
Vulkan implementation of ocean waves in both **DFT** and **FFT**

It is now basically a well organized version of Discrete Fourier Transformation Ocean Wave. 
The base project (base.lib) from the SaschaWillems examples codebase was used as a framework to provide a wrapper of Vulkan APIs.

You can use mouse and wheel to navigate. For now I do not set any cull mode.

#### Performance
It is very slowly due to the algorithm - Discrete Fourier Transform **(DFT)**.
In the latest version I update it to the Fast Fourier Transform algorithm **(FFT)**.

|  Mode   |   FPS   |                    Device                     |
| :-----: | :-----: | :-------------------------------------------: |
|  Debug  |  120+   | GPU(RTX 3050Ti Laptop) / CPU(Intel i7-11800H) |
| Release |  300~   | GPU(RTX 3050Ti Laptop) / CPU(Intel i7-11800H) |
|  Debug  |  120+   |   GPU(RTX 3070 Laptop) / CPU(AMD R7 5800H)    |
| Release | 680~950 |   GPU(RTX 3070 Laptop) / CPU(AMD R7 5800H)    |

#### Notice 
You can select your GPU in the line 870 in ***.\base01\vulkanexamplebase.cpp***. 
```
869     ...
870     uint32_t selectedDevice  = 0;
871     ...
```
The variable ***selectedDevice*** is 1 for default. You can change it to 0 in case:
- `You want to switch to another GPU` 
- `You only have only one GPU and the program crashes at the GPU selection.`
 
![wave](https://user-images.githubusercontent.com/42757470/159940977-d4095f51-8ebb-47dd-9f47-5e612a761801.png)
~ (The GPU will show on the window title and the control panel as well.)
