


## 使い方

```python
from percival.additional_calculation.domain.additional_calculation import AdditionalCalculation, AdditionalCalculationDirectory
from percival.additional_calculation.domain.value_objects import CalculationMethod, GaussianLog


log = GaussianLog.parse_file("tests/resources/entries/entry1/entry1_confs_0.out")
method = CalculationMethod("p tddft")
additional_calc = AdditionalCalculation(log, method)
acdict = AdditionalCalculationDirectory("add_calc", additional_calc)
acdict.additional_calculation.generate_gaussian_input(acdict.directory_path + "/fff.gjf")

```


## 概要

**Gaussianで計算したログ[1]** を指定して、 **追加の計算[2]** を行う。

**追加の計算[2]** は元の **計算ログ[1]** と同じ階層の **追加計算用ディレクトリ[2]** を作ってその中で行う。

**追加の計算[2]** は **計算メソッド[4]** を指定して、追加計算用インプットファイルを生成する。


### ユビキタス言語単語帳

1. GaussianLog
    - 量子化学計算のデファクトスタンダードになってるソフトウェアのログ
2. AdditionalCalculation
    - 構造最適化などを行なったのち、多くの場合、その構造を元に追加の計算をする。
    - このコンテキストではこの追加計算のインプットを適切に計算することを目的にする。
3. AdditionalCalculationDirectory
    - 元となる構造と同じ階層のディレクトリを作って、その中に追加計算のインプットや結果を入れておくと分かりやすい。
4. CalculationMethod
    - 追加計算用のメソッド
    - 多くの場合精度をあげたり、溶媒効果を取り入れた計算を行う。