# Vericrodo JPBC project

本项目依据 `Vericrodo_Main_Document(5).tex` 实现方案的前两个算法，并沿用
`CDECon.zip` 中 `InitGlob.java`、`InitAuth.java` 的 JPBC 编程和性能统计形式。

## 已实现算法

1. `InitGlob.java`：实现
   `InitGlob(kappa, Theta, U) -> pp`，其中
   `pp = {p, e, g1, g2, g3, F(.), H(.)}`。
2. `InitAuth.java`：实现
   `InitAuth(pp, sigma) -> (Apk_sigma, Ask_sigma)`，其中
   `Apk_sigma = {e(g1,g2)^alpha_sigma, g2^x_sigma}`，
   `Ask_sigma = {g1^alpha_sigma, x_sigma}`。

`InitAuth` 还包含公开密钥与秘密密钥的一致性检查：

- `e(g1^alpha_sigma, g2) = e(g1,g2)^alpha_sigma`；
- `g2^x_sigma` 等于公开密钥中的第二个分量。

## 曲线配置

两个算法默认统一加载 `d159.properties`，避免公共参数生成与权威初始化使用
不同曲线。项目同时原样保留了 CDECon 工程中的全部曲线配置文件，运行时可通过
第一个命令行参数切换，例如 `g149.properties`。

## JPBC 依赖

`Vericrodo.iml` 沿用 CDECon 的依赖目录结构，默认从项目同级目录读取：

```text
../Ideal_JPBC_JDK/jpbc-2.0.0/jars/
```

在 IntelliJ IDEA 中打开项目后，若本机 JPBC 目录不同，只需在 Project
Structure 中重新指定对应的 JPBC 2.0.0 JAR 文件。

## 运行参数

两个类均可直接运行，命令行参数格式相同：

```text
args[0] = 曲线配置文件路径，默认 d159.properties
args[1] = 重复次数，默认 100
args[2] = kappa，默认 128
```

`InitGlob` 按 `|U| = 10,20,...,100` 和 `|Theta| = 3,6,...,30` 输出平均计算
时间与公共参数存储开销；`InitAuth` 按 `|Theta| = 3,6,...,30` 输出全部权威初始化
的平均计算时间及 `Apk`、`Ask` 存储开销。
