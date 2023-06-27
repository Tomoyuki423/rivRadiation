# InfernoMelter 仕様書

Tomoyuki Kawashima

## この文章の概要

## InfernoMelter概要

### 開発目的

## ファイル構成

## コンパイルについて

## 名前空間

### 5.1 概要

本節では、名前空間の構造について述べる。FDPS は ParticleSimulator という名前空間で 囲まれている。以下では、ParticleSimulator 直下にある機能と、ParticleSimulator にネスト されている名前空間について述べる。 

### 5.2 ParticleSimulator

FDPS の標準機能すべては名前空間 ParticleSimulator の直下にある。 名前空間 ParticleSimulator は以下のように省略されており、この文書におけるあとの記述 でもこの省略形を採用する。 ✓ ✏ namespace PS = ParticleSimulator; ✒ ✑ 名前空間 ParticleSimulator の下にはいくつかの名前空間が拡張機能毎にネストされてい る。拡張機能には ParticleMesh がある。以下では拡張機能の名前空間について記述する

## データ型

### 概要

FDPSでは独自の整数型、実数型、ベクトル型、対称行列型、PS::SEARCH MODE 型、列 挙型が定義されている。整数型、実数型、ベクトル型、対称行列型に関しては必ずしもここに 挙げるものを用いる必要はないが、これらを用いることを推奨する。PS::SEARCH MODE 型、列挙型は必ず用いる必要がある。以下、整数型、実数型、ベクトル型、対称行列型、 PS::SEARCH MODE 型、列挙型の順に記述する。

### 6.4.3 IM::Vector3

PS::Vecotr3 は x, y, z の 3 要素を持つ。これらに対する様々な API や演算子を定義した。 それらの宣言を以下に記述する。この節ではこれらについて詳しく記述する。 ソースコード 9: Vector3 1 namespace ParticleSimulator { 2 template < typename T > 3 class Vector3 { 4 public : 5 //メンバ変数は以下の二つのみ。 6 T x , y , z ; 7 8 //コンストラクタ 9 Vector3 () : x ( T (0)) , y ( T (0)) , z ( T (0)) {} 10 Vector3 ( const T _x , const T _y , const T _z ) : x ( _x ) , y ( _y ) , z( _z ) {} 11 Vector3 ( const T s ) : x (s ) , y ( s ), z ( s ) {} 12 Vector3 ( const Vector3 & src ) : x ( src . x ) , y ( src . y ) , z ( src .z ) {} 13 14 //代入演算子 15 const Vector3 & operator = ( const Vector3 & rhs ); 16 17 // [] 演 算 子 18 const T & opertor []( const int i ); 19 T & operator []( const int i ); 20 21 //加減算 22 Vector3 operator + ( const Vector3 & rhs ) const ; 23 const Vector3 & operator += ( const Vector3 & rhs ); 24 Vector3 operator - ( const Vector3 & rhs ) const ; 25 const Vector3 & operator -= ( const Vector3 & rhs ); 26 27 //ベクトルスカラ積 28 Vector3 operator * ( const T s ) const ; 29 const Vector3 & operator *= ( const T s ); 30 friend Vector3 operator * ( const T s , const Vector3 & v ); 31 Vector3 operator / ( const T s ) const ; 32 const Vector3 & operator /= ( const T s ); 37 33 34 //内積 35 T operator * ( const Vector3 & rhs ) const ; 36 37 //外積( 返 り 値 は ス カ ラ! ! ) 38 T operator ^ ( const Vector3 & rhs ) const ; 39 40 // Vector3  へ の 型 変 換 41 template < typename U > 42 operator Vector3  () const ; 43 }; 44 } 6.4.3.1 コンストラクタ ✓ ✏ template PS::Vector3::Vector3() ✒ ✑ • 引数 なし。 • 機能 デフォルトコンストラクタ。メンバ x,y は 0 で初期化される。 ✓ ✏ template PS::Vector3::Vector3(const T _x, const T _y) ✒ ✑ • 引数 x: 入力。const T 型。 y: 入力。const T 型。 • 機能 メンバ x、y をそれぞれ x、 y で初期化する。 ✓ ✏ template PS::Vector3::Vector3(const T s); ✒ ✑ 38 • 引数 s: 入力。const T 型。 • 機能 メンバ x、y を両方とも s の値で初期化する。 6.4.3.2 コピーコンストラクタ ✓ ✏ template PS::Vector3::Vector3(const PS::Vector3 & src) ✒ ✑ • 引数 src: 入力。const PS::Vector3 &型。 • 機能 コピーコンストラクタ。src で初期化する。 6.4.3.3 メンバ変数 ✓ ✏ template T PS::Vector3::x; template T PS::Vector3::y; template T PS::Vector3::z; ✒ ✑ • 機能 メンバ x、y、z を直接操作出来る。 6.4.3.4 代入演算子 ✓ ✏ template const PS::Vector3 & PS::Vector3::operator = (const PS::Vector3 & rhs); ✒ ✑ 39 • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 const PS::Vector3 &型。rhs の x,y の値を自身のメンバ x,y に代入し自身の参照を 返す。代入演算子。 6.4.3.5 [] 演算子 ✓ ✏ template const T & PS::Vector3::operator[] (const int i); ✒ ✑ • 引数 i: 入力。const int 型。 • 返り値 const  &型。ベクトルの i 成分を返す。 • 備考 直接メンバ変数を指定する場合に比べ、処理が遅くなることがある。 ✓ ✏ template T & PS::Vector3::operator[] (const int i); ✒ ✑ • 引数 i: 入力。const int 型。 • 返り値 &型。ベクトルの i 成分を返す。 • 備考 直接メンバ変数を指定する場合に比べ、処理が遅くなることがある。 40 6.4.3.6 加減算 ✓ ✏ template PS::Vector3 PS::Vector3::operator + (const PS::Vector3 & rhs) const; ✒ ✑ • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 PS::Vector3 型。rhs の x,y の値と自身のメンバ x,y の値の和を取った値を返す。 ✓ ✏ template const PS::Vector3 & PS::Vector3::operator += (const PS::Vector3 & rhs); ✒ ✑ • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 const PS::Vector3 &型。rhs の x,y の値を自身のメンバ x,y に足し、自身を返す。 ✓ ✏ template PS::Vector3 PS::Vector3::operator - (const PS::Vector3 & rhs) const; ✒ ✑ • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 PS::Vector3 型。rhs の x,y の値と自身のメンバ x,y の値の差を取った値を返す。 ✓ ✏ template const PS::Vector3 & PS::Vector3::operator -= (const PS::Vector3 & rhs); ✒ ✑ 41 • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 const PS::Vector3 &型。自身のメンバ x,y から rhs の x,y を引き自身を返す。 6.4.3.7 ベクトルスカラ積 ✓ ✏ template PS::Vector3 PS::Vector3::operator * (const T s) const; ✒ ✑ • 引数 s: 入力。const T 型。 • 返り値 PS::Vector3 型。自身のメンバ x,y それぞれに s をかけた値を返す。 ✓ ✏ template const PS::Vector3 & PS::Vector3::operator *= (const T s); ✒ ✑ • 引数 rhs: 入力。const T 型。 • 返り値 const PS::Vector3 &型。自身のメンバ x,y それぞれに s をかけ自身を返す。 ✓ ✏ template PS::Vector3 PS::Vector3::operator / (const T s) const; ✒ ✑ • 引数 s: 入力。const T 型。 • 返り値 PS::Vector3 型。自身のメンバ x,y それぞれを s で割った値を返す。 42 ✓ ✏ template const PS::Vector3 & PS::Vector3::operator /= (const T s); ✒ ✑ • 引数 rhs: 入力。const T 型。 • 返り値 const PS::Vector3 &型。自身のメンバ x,y それぞれを s で割り自身を返す。 6.4.3.8 内積、外積 ✓ ✏ template T PS::Vector3::operator * (const PS::Vector3 & rhs) const; ✒ ✑ • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 T 型。自身と rhs の内積を取った値を返す。 ✓ ✏ template T PS::Vector3::operator ^ (const PS::Vector3 & rhs) const; ✒ ✑ • 引数 rhs: 入力。const PS::Vector3 &型。 • 返り値 T 型。自身と rhs の外積を取った値を返す。 6.4.3.9 Vector3 への型変換 ✓ ✏ template template  PS::Vector3::operator PS::Vector3 () const; ✒ ✑ 43 • 引数 なし • 返り値 const PS::Vector3 型。 • 機能 const PS::Vector3 型を const PS::Vector3 型にキャストする。