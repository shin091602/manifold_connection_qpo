# manifold_connection_qpo リポジトリ説明

## 概要
このリポジトリは、軌道力学における**制限三体問題（Circular Restricted Three-Body Problem: CR3BP）**の数値計算とマニフォールド解析を行うMATLABコードです。特に**準周期軌道（Quasi-Periodic Orbit: QPO）**と**準周期不変トーラス（Quasi-Periodic invariant Tori: QPT）**の計算・可視化・解析に特化した研究用ソフトウェアです。

## 主要機能
### 1. 軌道計算
- 制限三体問題の軌道積分
- 周期軌道の計算と修正
- 準周期軌道の継続計算
- ハミルトン系の数値解法

### 2. マニフォールド解析
- 不変マニフォールドの計算
- トーラス構造の数値表現
- GMOSアルゴリズムによる準周期不変トーラスの計算
- マニフォールド接続問題の解析

### 3. 可視化・プロット機能
- 軌道の3D可視化
- 不変トーラスのプロット
- モノドロミー行列の固有値解析
- フーリエ解析による軌道特性の可視化

## ディレクトリ構成

### `/Functions/`
基本的な数値計算関数群
- `fun_cr3bp.m`: 制限三体問題の運動方程式
- `parameter.m`: 天体系パラメータ（Sun-Earth, Earth-Moon等）
- `librationPoints.m`: ラグランジュ点の計算
- `fun_manifold_cr3bp.m`: マニフォールド計算関数
- `fun_stm_cr3bp.m`: 状態遷移行列の計算

### `/function_QPT/`
準周期軌道・トーラス計算の専用関数
- `/CR3BP/`: 制限三体問題用の特殊化された関数群
- `F_qpoms.m`, `DF_qpoms.m`: QPOMS（準周期軌道多重射撃法）関数
- `PAC_qpoms.m`: PAC（擬似弧長継続法）実装
- `fourier_matrix.m`: フーリエ変換マトリックス

### `/Ayano functions/`
Ayano Tsuruta氏による拡張関数群
- `/A/`: 行列計算関数（2次元・4次元対応）
- `/ODE/`: 常微分方程式ソルバー（自然・非自然座標系）
- `/plot/`: 専用プロット関数

### `/susumu functions/`
susumu氏による補間・解析関数
- `direction_manifold_interpolation.m`: マニフォールド方向の補間
- `fun_Fourier_interpolation.m`: フーリエ補間関数

## メインスクリプト

### `QPO_gmos_CR3BP_original.m`
GMOSアルゴリズムを用いた2次元準周期不変トーラス計算のメインスクリプト
- Nicola Baresi et al.の論文に基づく実装
- 完全数値的な準周期不変トーラス族の継続計算

### `QPO_gmos_CR3BP_2D_matrix.m`
行列形式での2次元トーラス計算の最適化版

### `hetero_fmincon_plane.m`
異質接続問題の最適化計算

## 対応天体系
- Sun-Earth系
- Earth-Moon系  
- Mars-Phobos系
- Mars-Deimos系

## 数値手法
- **GMOS (Galerkin Method on Orbit Segments)**: 軌道セグメントに対するガラーキン法
- **PAC (Pseudo-Arc-length Continuation)**: 擬似弧長継続法
- **Multiple Shooting**: 多重射撃法
- **Differential Correction**: 微分修正法

## 参考文献
- Nicola Baresi et al., "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodynamics"
- Soi Yamaguchi (2023年度修士論文)
- Damennick Bolte Henry (Julia実装版を参考)

## 使用方法
1. MATLABで該当ディレクトリを開く
2. `parameter.m`で計算対象の天体系を設定
3. メインスクリプト（`QPO_gmos_*.m`）を実行
4. 結果は構造体`sol_qpos`に格納される

## 注意事項
- MATLAB R2020a以降を推奨
- Symbolic Math Toolboxが必要な関数があります
- 大規模計算では十分なメモリ容量を確保してください
- `.gitignore`により、`.mat`ファイルや`.fig`ファイルは管理対象外です

## 開発・貢献
ブランチ命名規則：`purpose/descriptive-name`形式
- `feature/`: 新機能追加
- `bugfix/`: バグ修正
- `enhancement/`: 精度・性能向上
- `experiment/`: 検証・実験用

詳細は`README.md`のBranch Naming Conventionセクションを参照してください。