# plot_manifold_qpo

for plotting manifold from a point on any surface of qpo.

### `.gitignore` の説明

このリポジトリでは、MATLAB による研究・開発において不要なファイルをバージョン管理から除外するために `.gitignore` を設定しています。

#### 除外対象の主なファイル・ディレクトリとその理由

- `*.asv`, `*.m~`, `*.autosave`  
  → MATLAB の自動保存ファイルや一時ファイル。ユーザーごとに生成されるため、管理の対象外とします。

- `*.mat`, `*.fig`  
  → 実験データや可視化ファイルなどの中間生成物。再現性が担保される限り、リポジトリでは管理しません。

- `slprj/`, `sfprj/`, `*_grt_rtw/`, `*_ert_rtw/`, `*.slxc`  
  → Simulink によるコード生成やシミュレーション結果。ビルドごとに変更されるため、除外しています。

- `.DS_Store`, `Thumbs.db`, `*.swp` など  
  → OS やエディタ特有のキャッシュ・設定ファイルで、プロジェクト本体には無関係です。

#### 注意事項

`.gitignore` に含まれているファイルの中でも、**重要な研究成果や共有すべきデータ**がある場合は、個別に明示し、必要に応じて `README.md` やドキュメントで補足してください。

## 🌿 Branch Naming Convention (目的ベース)

本リポジトリでは、作業の**目的に基づいてブランチを命名**します。命名ルールを以下に示します。

### 📌 命名フォーマット

- `purpose`：作業の目的（カテゴリ）
- `descriptive-name`：作業内容の要約（ケバブケースを推奨）

---

### 🔧 Purpose 一覧

| Purpose       | 説明                             | 例                                        |
| ------------- | -------------------------------- | ----------------------------------------- |
| `feature`     | 新機能追加                       | `feature/interpolated-point-management`   |
| `bugfix`      | バグ修正                         | `bugfix/nullspace-tolerance-adjustment`   |
| `refactor`    | コード構造の改善（機能変更なし） | `refactor/pac-function-modularization`    |
| `enhancement` | 精度改善・性能向上               | `enhancement/improve-torus-plot-accuracy` |
| `experiment`  | 検証・実験用のブランチ           | `experiment/test-new-integration-method`  |
| `doc`         | ドキュメント追加・修正           | `doc/add-jacobi-constant-description`     |
| `test`        | テストコード追加・修正           | `test/interpolation-function-unit-tests`  |
| `hotfix`      | 緊急修正（公開後の重大バグ等）   | `hotfix/fix-stm-crash`                    |

---

### ✅ 命名例

- `feature/qpt-continuation`
- `bugfix/fix-eigenvalue-sorting`
- `refactor/clean-gmos-interface`
- `enhancement/tune-svd-threshold`
- `doc/update-usage-section`
- `test/add-po-validation-tests`

---

### 💡 Tips

- **小文字・ケバブケース（`-`）** を推奨
- `descriptive-name` はできるだけ**簡潔かつ具体的に**
- Issue に関連づける場合は末尾に `-#番号` をつけると便利：
  - 例：`bugfix/fix-nullspace-eval-#42`

---
