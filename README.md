# cpp-tsp3
使用遗传算法解决 tsp 问题，现代软件工程大作业

- 项目使用 c++ 编写，附带使用 javascript 编写的仿真程序
- c++使用了面向对象的方式来组织，附带注释，可读性应该会好一点
- 使用了一些常量来控制计算过程，旁边都有注释
- 图方便就把c++代码写在一个文件里了
- 项目使用 vs2017 编译，建议在运行的时候切换到 Release 模式，速度和 Debug 模式不在一个数量级上
- 也可以使用 g++ 直接编译 `cpp-tsp3.cpp` 文件，记得把那句 `#include "stdafx.h"` 去掉
- 使用了 c++ 11 的语法和库，请确认编译器版本
- 具体的使用文档在 docx 文件里

另外，这次尝试了完全不使用指针或手动在堆上分配内存来编写代码，使用 move 语义应该能让效率和使用指针的场景差不多。
