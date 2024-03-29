# Japanese
## v4_nays3dv
【あっと驚くスーパー3次元密度流計算モデル】
Nays3Dvは自由水面を持つ3次元流れ場とその中の密度変化を伴う流れ場計算するソルバである。3次元流れの運動方程式、連続式および物質濃度の移流拡散方程式は非定常、境界適合の条件で計算される。
2020年から2021年の世界的なコロナ禍の中、普段であれば世界中あちこち遊び回っている札幌のS氏が、緊急事態宣言などの影響でほぼ自宅軟禁状態にされた状態で開発されたモデルである。
密度流は、流体中の温度、物質濃度などの違いによって発生する浮力、重力流等の効果が反映出来るようになっている。
河川河口付近の塩水遡上、ダム湖内の温度・密度流循環などに使用可能である。
また、密度流の特性を3次元的にヴィジュアルに表現出来るので、教育用ツールとしても良いかも知れない。

事例集はこちらを参照ください。
https://i-ric.org/yasu/Nays3dv_JP/index.html

## リリースノート
### ver.1.7.24012401
* 計算中の計算結果の読み込み機能に軽微な修正が行われました。

### ver.1.7.24011001
* 開発中のため使用できない機能の選択が可能だった問題を修正しました。

### ver.1.7.23021302
* READMEに事例集のリンクを追加。
### ver.1.7.23021301
* ライセンスの情報を修正
### ver.1.7.22101301
* iRIC v4版をリリース

# English

## v4_nays3dv
Nays3Dv is 3Dimensional model developed for calculation of vertical and horizontal movement of fluid with density currents.
The Nays3Dv solver is developed by Professor Yasuyuki Shimizu from Hokkaido University.
Density currents occured due to density differences, arise from temperature variations, suspended solids or dissolved materials. Therefore, formation and evolution of density currents are induced by natural conditions such as saline intrusions, oil spills etc. Density flow is important for problems in lakes, reservoirs and estuaries.

Check here for case studies.
https://i-ric.org/yasu/Nays3dv/index.html

## Release notes
ver.1.7.24012401
* Minor adjustments to the result loading function during calculations.
### ver.1.7.24011001
* Addressed the issue where it was possible to select features that are not available due to being in development.
### ver.1.7.23021302
* Add case study link to README.
### ver.1.7.23021301
* Fix license information
### ver.1.7.22101301
* Released for iRIC v4