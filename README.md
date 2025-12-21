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

## ⚠️ 已知问题

**支路无功功率计算误差**：与 MATPOWER 标准结果对比，程序输出的支路无功功率（Qij 和 Qji）存在较大误差，误差范围约在 12%-28% 之间。该误差可能源于支路功率计算中对并联电纳（线路充电电纳）的处理方式或计算顺序的差异。

**建议**：
- ✅ 节点电压、有功功率和损耗计算结果准确，可放心使用
- ⚠️ 支路无功功率结果仅供参考，使用时需谨慎验证
- 🔧 如需高精度支路无功功率，建议参考 MATPOWER 或其他成熟潮流计算软件

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
│   ├── case4.m
│   ├── case5.m
│   ├── case30.m
│   ├── case118.m
│   └── case2383.m
└── README.md                  # 项目说明文档
```

## 🙏 致谢

本项目在开发过程中参考了多个开源潮流计算项目的实现思路，特别感谢相关开发者的贡献。希望这份代码能够帮助到正在学习电力系统潮流计算的同学们。

## 📄 许可证

本项目采用 MIT 许可证，详见 [LICENSE](LICENSE) 文件。

---

## ⚠️ 重要提示

1. **本项目为课程设计作业，仅供学习交流使用**。如有问题或建议，欢迎提出 Issue 或 Pull Request。

2. **支路无功功率计算存在误差**：程序输出的支路无功功率（Qij 和 Qji）与 MATPOWER 标准结果存在较大误差（约 12%-28%），使用时请谨慎参考。节点电压、有功功率和损耗计算结果准确可靠。
