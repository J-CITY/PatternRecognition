package Lab2No1

import (
	"fmt"
	"math"
	"math/rand"
)

//lab2No1
type Lab2No1 struct {
	N int
	X []int

	Path string

	Frequencies []int
	D float64
	Sigma float64
	MiddleX float64
	Asymmetry float64
	Kurtosis float64

	FXmin float64
	FXmax float64
	HXmin float64
	HXmax float64

}
func (in *Lab2No1) Init() {
	in.N = 0
	in.FXmin = 0.0
	in.FXmax = 0.0
	in.HXmin = 0.0
	in.HXmax = 0.0
}
func (in *Lab2No1) CalcParams() {
	for i := 0; i < in.N; i+=1 {
		in.MiddleX += float64(in.X[i])
	}
	in.MiddleX /= float64(in.N)

	for i := 0; i < in.N; i+=1 {
		in.D += (float64(in.X[i])-in.MiddleX)*(float64(in.X[i])-in.MiddleX)
	}
	in.D /= float64(in.N)
	in.Sigma = math.Sqrt(in.D)

	sumX := 0.0
	for i := 0; i < in.N; i+=1 {
		tmp := float64(in.X[i])-in.MiddleX
		sumX +=  math.Pow(tmp, 3)
	}
	mu3 := sumX / float64(in.N)
	in.Asymmetry = mu3 / math.Pow(math.Sqrt(in.D), 3)

	sumX = 0.0
	for i := 0; i < in.N; i+=1 {
		tmp := float64(in.X[i])-in.MiddleX
		sumX += math.Pow(tmp, 4)
	}
	mu4 := sumX / float64(in.N)
	in.Kurtosis = mu4 / math.Pow(math.Sqrt(in.D), 4) - 3

	fmt.Println("Params:")
	fmt.Println(in.MiddleX)
	fmt.Println(in.D)
	fmt.Println(in.Asymmetry)
	fmt.Println(in.Kurtosis)
}
func (in *Lab2No1) GenX() {
	in.X = nil
	//for i := 0; i < len(in.Frequencies); i+=1 {
	//	in.Frequencies[i] = 0
	//}
	in.Frequencies = nil
	in.Frequencies = make([]int, 100)
    for i := 0; i < in.N; i+=1 {
		j := 0
		event := 0
		for j < 1 {
			v := rand.Intn(4)
			if v == 0 {
				in.X = append(in.X, event)
				break
			} else {
				event += 1
			}
		}
		in.Frequencies[event] += 1
	}
	in.CalcParams()
	for i := 0; i < len(in.Frequencies); i+=1 {
		//in.Frequencies[i] /= in.N
		fmt.Println("F[", i+1, "]: ", in.Frequencies[i])
	}
}

func (in *Lab2No1) SetX(X []byte) {
	in.N = len(X)
	in.X = nil
	for i := 0; i < in.N; i+=1 {
		in.X = append(in.X, int(X[i]))
	}

	//for i := 0; i < len(in.Frequencies); i+=1 {
	//	in.Frequencies[i] = 0
	//}
	
	in.Frequencies = nil
	in.Frequencies = make([]int, 100)
    for i := 0; i < in.N; i+=1 {
		in.Frequencies[X[i]] += 1
	}
	in.CalcParams()
	for i := 0; i < len(in.Frequencies); i+=1 {
		fmt.Println("F[", i+1, "]: ", in.Frequencies[i])
	}
}
