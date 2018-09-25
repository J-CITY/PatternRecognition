package Lab2No2

import (
	"fmt"
	"math"
	"math/rand"
)

//lab2No2
type Lab2No2 struct {
	N int
	X []float64

	Path string

	Frequencies []int
	FN int

	MiddleX float64
	D float64
	Sigma float64
	Asymmetry float64
	Kurtosis float64
	
	a float64
	k float64
	x0, xN, c float64

	FXmin float64
	FXmax float64
	HXmin float64
	HXmax float64

}
func (in *Lab2No2) Init() {
	in.N = 0
	in.a = 0.5
	in.k = 1
	in.x0 = 0
	in.xN = 1
	in.c = 100
	in.FXmin = 0.0
	in.FXmax = 0.0
	in.HXmin = 0.0
	in.HXmax = 0.0
}

func (in *Lab2No2) SetX0(_x float64) {
	in.x0 = _x
}
func (in *Lab2No2) SetXN(_x float64) {
	in.xN = _x
}
func (in *Lab2No2) SetC(_x float64) {
	in.c = _x
}

func (in *Lab2No2) GetX0() float64 {
	return in.x0
}
func (in *Lab2No2) GetXN() float64 {
	return in.xN
}
func (in *Lab2No2) GetC() float64 {
	return in.c
}

func (in *Lab2No2) CalcParams() {
	for i := 0; i < in.N; i+=1 {
		in.MiddleX += in.X[i]
	}
	in.MiddleX /= float64(in.N)

	for i := 0; i < in.N; i+=1 {
		in.D += (in.X[i]-in.MiddleX)*(in.X[i]-in.MiddleX)
	}
	in.D /= float64(in.N)
	in.Sigma = math.Sqrt(in.D)

	sumX := 0.0
	for i := 0; i < in.N; i+=1 {
		tmp := in.X[i]-in.MiddleX
		sumX +=  math.Pow(tmp, 3)
	}
	mu3 := sumX / float64(in.N)
	in.Asymmetry = mu3 / math.Pow(math.Sqrt(in.D), 3)

	sumX = 0.0
	for i := 0; i < in.N; i+=1 {
		tmp := in.X[i]-in.MiddleX
		sumX += math.Pow(tmp, 4)
	}
	mu4 := sumX / float64(in.N)
	in.Kurtosis = mu4 / math.Pow(math.Sqrt(in.D), 4) - 3

	fmt.Println(in.MiddleX)
	fmt.Println(in.D)
	fmt.Println(in.Asymmetry)
	fmt.Println(in.Kurtosis)
}


func (in *Lab2No2)getFX(x float64) float64 {
    return in.a*math.Exp(-in.k*x)*x*x;
}
func (in *Lab2No2) GenX(_N int) {
	customRand := func() int {
		_a, _b, _c := in.x0, in.xN, in.c // 0, 1, 1 for default
    	r1 := rand.Float64()
    	r2 := rand.Float64()

    	x1 := _a + (_b-_a)*r1;
    	y1 := _c*r2;
    	if y1 < in.getFX(x1) {
			in.X = append(in.X, x1)
    	    return 0
		} else {
    	    return -1
    	}
	}
	in.X = nil
	i := 0
    for i < in.N {
        if customRand() == 0 {
            i+=1
        }
	}
	if _N > 0 {
		in.FN = _N
	} else {
		in.FN = int(5*math.Log(float64(in.N)))
		//fmt.Println(in.FN)
	}
	//in.FN = math.sqrt(in.N)
	//in.FN = 4*math.log(in.N)
	//in.FN = 5*math.log(float64(in.N)/10)
	step := float64(in.xN-in.x0)/float64(in.FN)
	
	in.Frequencies = nil
	in.Frequencies = make([]int, in.FN)

	for i := 0; i < in.N; i+=1 {
		j := 0
		for in.X[i] > in.x0 + float64(j)*step {
			j+=1
		}
		j-=1
		//fmt.Println(j)
		in.Frequencies[j]+=1
	}
	in.CalcParams()
	for i := 0; i < len(in.Frequencies); i+=1 {
		fmt.Println("F[", i+1, "]: ", in.Frequencies[i])
	}
}

func (in *Lab2No2) SetX(X []float64, _N int) {
	in.N = len(X)
	in.X = nil
	for i := 0; i < in.N; i+=1 {
		in.X = append(in.X, X[i])
	}

	
	if _N > 0 {
		in.FN = _N
	} else {
		in.FN = int(5*math.Log(float64(in.N)))
	}
	in.Frequencies = nil
	in.Frequencies = make([]int, in.FN)
	//in.FN = math.sqrt(in.N)
	//in.FN = 4*math.log(in.N)
	//in.FN = 5*math.log(float64(in.N)/10)
	step := float64(in.xN-in.x0)/float64(in.FN)
	for i := 0; i < in.N; i+=1 {
		j := 0
		for in.X[i] > in.x0 + float64(j)*step {
			j+=1
		}
		j-=1
		in.Frequencies[j]+=1
	}
	in.CalcParams()
	for i := 0; i < len(in.Frequencies); i+=1 {
		fmt.Println("F[", i+1, "]: ", in.Frequencies[i])
	}
}
