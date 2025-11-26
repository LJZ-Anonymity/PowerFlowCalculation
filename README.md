# 电力系统潮流计算

> 基于牛顿-拉夫逊法的电力系统稳态潮流计算程序（课程设计）

[![最后更新](https://img.shields.io/github/last-commit/LJZ-Anonymity/PowerFlowCalculation)](https://github.com/LJZ-Anonymity/PowerFlowCalculation/commits)
[![开源协议](https://img.shields.io/badge/License-GNU%20AGPL%20v3.0-blue.svg)](https://github.com/LJZ-Anonymity/PowerFlowCalculation/blob/main/LICENSE)
[![代码行数](https://aschey.tech/tokei/github/LJZ-Anonymity/PowerFlowCalculation)](https://github.com/LJZ-Anonymity/PowerFlowCalculation)
[![主要语言](https://img.shields.io/github/languages/top/LJZ-Anonymity/PowerFlowCalculation)](https://github.com/LJZ-Anonymity/PowerFlowCalculation)

## 📋 项目简介

本项目是一个使用 C 语言实现的电力系统潮流计算程序，采用**直角坐标牛顿-拉夫逊法**求解电力系统稳态潮流问题。程序支持 PQ 节点、PV 节点和平衡节点的处理，能够计算从简单四节点系统到大型电力系统的潮流分布。

## ✨ 功能特性

- ✅ **直角坐标牛顿-拉夫逊法**：使用 e、f 变量进行迭代求解
- ✅ **完整的节点类型支持**：PQ 节点、PV 节点、平衡节点
- ✅ **节点导纳矩阵自动形成**：从线路参数自动构建 G+jB 矩阵
- ✅ **高精度计时**：使用 Windows 高精度计时器测量计算时间
- ✅ **标准格式输出**：按照 MATPOWER 的 `printpf` 风格输出结果
- ✅ **详细的计算报告**：包含系统摘要、节点数据、支路数据及总计信息

## 🔧 技术实现

- **编程语言**：C 语言
- **开发环境**：Windows（Visual Studio）
- **算法**：牛顿-拉夫逊迭代法（直角坐标形式）
- **线性求解**：高斯消元法（带列主元选择）
- **矩阵运算**：自主实现，不依赖外部矩阵库

## 📁 文件结构

```
PowerFlowCalculation/
├── PowerFlowCalculation.c    # 主程序源代码
├── Input/                     # 输入数据目录
│   ├── case4.txt
│   ├── case5.txt
│   ├── case30.txt
│   ├── case118.txt
│   └── case2383.txt
└── README.md                  # 项目说明文档
```

## 📊 支持的系统规模

- ✅ case4（4 节点系统）
- ✅ case5（5 节点系统）
- ✅ case30（30 节点系统）
- ✅ case118（118 节点系统）
- ✅ case2383（2383 节点系统）

## 🙏 致谢

本项目在开发过程中参考了多个开源潮流计算项目的实现思路，特别感谢相关开发者的贡献。希望这份代码能够帮助到正在学习电力系统潮流计算的同学们。

## 📄 许可证

本项目采用 GNU AGPL v3.0 许可证，详见 [LICENSE](LICENSE) 文件。

---

**注意**：本项目为课程设计作业，仅供学习交流使用。如有问题或建议，欢迎提出 Issue 或 Pull Request。
