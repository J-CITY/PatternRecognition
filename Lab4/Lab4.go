package Lab4

import (
	//"fmt"
	"github.com/mjibson/go-dsp/fft"
)
import "github.com/mjibson/go-dsp/spectral"
import (
	"image/color"
	//"math/rand"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	//"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"math"
	"time"
	"math/rand"

	"math/cmplx"
	//"gonum.org/v1/plot/vg/draw"
)
	
type Signal interface {
    GetMax() float64
    GetMin() float64
	
	
	Run(in float64) float64
}


type DelayedSingleImpulse struct {
    delay float64
}
func (s *DelayedSingleImpulse) SetDelay(in float64) {
    s.delay = in
}
func (s *DelayedSingleImpulse) SetParams(in float64) {
    s.delay = in
}
func (s *DelayedSingleImpulse) GetMax() float64 {
    return 1
}
func (s *DelayedSingleImpulse) GetMin() float64 {
    return 0
}
func (s *DelayedSingleImpulse) Run(in float64) float64 {
	if in == s.delay {
		return 1
	}
	return 0
}

type DelayedUnitJump struct {
    delay float64
}
func (s *DelayedUnitJump) SetDelay(in float64) {
    s.delay = in
}
func (s *DelayedUnitJump) SetParams(in float64) {
    s.delay = in
}
func (s *DelayedUnitJump) GetMax() float64 {
    return 1
}
func (s *DelayedUnitJump) GetMin() float64 {
    return 0
}
func (s *DelayedUnitJump) Run(in float64) float64 {
	if in >= s.delay {
		return 1
	}
	return 0
}

type DiscretizedDecreasingExponent struct {
    a float64
}
func (s *DiscretizedDecreasingExponent) SetParams(in float64) {
    s.a = in
}
func (s *DiscretizedDecreasingExponent) SetA(in float64) {
    s.a = in
}
func (s *DiscretizedDecreasingExponent) GetMax() float64 {
    return 1
}
func (s *DiscretizedDecreasingExponent) GetMin() float64 {
    return 0
}
func (s *DiscretizedDecreasingExponent) Run(in float64) float64 {
	return math.Pow(s.a, in)
}

type SampledSine struct {
	amplitude float64
	frequency float64
    initialPhase float64
}
func (s *SampledSine) SetParams(inA, inF, inP float64) {
	s.amplitude = inA
	s.frequency = inF
	s.initialPhase = inP
}
func (s *SampledSine) SetAmplitude(in float64) {
    s.amplitude = in
}
func (s *SampledSine) SetFrequency(in float64) {
    s.frequency = in
}
func (s *SampledSine) SetInitialPhase(in float64) {
    s.initialPhase = in
}
func (s *SampledSine) GetMax() float64 {
    return s.amplitude
}
func (s *SampledSine) GetMin() float64 {
    return -s.amplitude
}
func (s *SampledSine) Run(in float64) float64 {
	return s.amplitude * math.Sin(in * s.frequency + s.initialPhase)
}

type Meander struct {
	period float64
}
func (s *Meander) SetParams(in float64) {
	s.period = in
}
func (s *Meander) SetPeriod(in float64) {
    s.period = in
}
func (s *Meander) GetMax() float64 {
    return 1
}
func (s *Meander) GetMin() float64 {
    return -1
}
func (s *Meander) Run(in float64) float64 {
	if math.Mod(in, s.period) < s.period / 2 {
		return 1
	}
	return -1
}

type Saw struct {
	period float64
}
func (s *Saw) SetParams(in float64) {
	s.period = in
}
func (s *Saw) SetPeriod(in float64) {
    s.period = in
}
func (s *Saw) GetMax() float64 {
    return 1
}
func (s *Saw) GetMin() float64 {
    return 0
}
func (s *Saw) Run(in float64) float64 {
	return  math.Mod(in, s.period) / s.period
}

type ExponentialEnvelope struct {
	amplitude float64
	frequency float64
	initialPhase float64
	envelopeWidth float64
}
func (s *ExponentialEnvelope) SetParams(a, f, i, e float64) {
	s.amplitude = a
	s.frequency = f
	s.initialPhase = i
	s.envelopeWidth = e
}
func (s *ExponentialEnvelope) GetMax() float64 {
    return 1
}
func (s *ExponentialEnvelope) GetMin() float64 {
    return 0
}
func (s *ExponentialEnvelope) Run(in float64) float64 {
	return  s.amplitude * math.Exp(-in / s.envelopeWidth) * math.Cos(in * s.frequency + s.initialPhase)
}

type BalancedEnvelope struct {
	amplitude float64
	frequency float64
	initialPhase float64
	envelopeFreq float64
}
func (s *BalancedEnvelope) SetParams(a, f, i, e float64) {
	s.amplitude = a
	s.frequency = f
	s.initialPhase = i
	s.envelopeFreq = e
}
func (s *BalancedEnvelope) GetMax() float64 {
    return 1
}
func (s *BalancedEnvelope) GetMin() float64 {
    return 0
}
func (s *BalancedEnvelope) Run(in float64) float64 {
	return s.amplitude * math.Cos(s.envelopeFreq*in) * math.Cos(s.frequency*in+s.initialPhase)
}

type TonalEnvelope struct {
	amplitude float64
	frequency float64
	initialPhase float64
	envelopeFreq float64
	depthIndex float64 //0..1
}
func (s *TonalEnvelope) SetParams(a, f, i, e, d float64) {
	s.amplitude = a
	s.frequency = f
	s.initialPhase = i
	s.envelopeFreq = e
	s.depthIndex = d
}
func (s *TonalEnvelope) GetMax() float64 {
    return 1
}
func (s *TonalEnvelope) GetMin() float64 {
    return 0
}
func (s *TonalEnvelope) Run(in float64) float64 {
	return s.amplitude * (1 + s.depthIndex * math.Cos(s.envelopeFreq*in)) * math.Cos(s.frequency*in + s.initialPhase)
}

var s1 = rand.NewSource(time.Now().UnixNano())
var r1 = rand.New(s1)

type WhiteNoise struct {
	IA, IB float64
}
func (s *WhiteNoise) SetParams(a, b float64) {
	s.IA = a
	s.IB = b
}
func (s *WhiteNoise) GetMax() float64 {
    return s.IB
}
func (s *WhiteNoise) GetMin() float64 {
    return s.IA
}
func (s *WhiteNoise) Run(in float64) float64 {
	return s.IA + r1.Float64()*(s.IB- s.IA)
}

type WhiteNoiseNorm struct {
	Med, D float64
	second float64
	ready bool

}
func (s *WhiteNoiseNorm) SetParams(m, d float64) {
	s.Med = m
	s.D = d

	s.ready = false;
	s.second = 0.0;
}
func (s *WhiteNoiseNorm) GetMax() float64 {
    return 1
}
func (s *WhiteNoiseNorm) GetMin() float64 {
    return 0
}

func (s *WhiteNoiseNorm) next() float64 {
	 if s.ready {
		s.ready = false;
		return s.second * s.D + s.Med
	 } else {
		u := 0.0
		v := 0.0
		ss := 0.0
		for (ss > 1.0 || ss == 0.0) {
			u = 2.0 * r1.Float64() - 1.0;
			v = 2.0 * r1.Float64() - 1.0;
			ss = u * u + v * v;
		}
		
		r := math.Sqrt(-2.0 * math.Log(ss) / ss);
		s.second = r * u;
		s.ready = true;
		return r * v * s.D + s.Med;
	 }
 }

func (s *WhiteNoiseNorm) Run(in float64) float64 {
	return s.next()
}

type Autoregr struct {
	A, B []float64

	hy []float64
	hx []float64
	second float64
	ready bool
	Med, D float64
}
func (s *Autoregr) SetParams(a, b []float64) {
	s.A = a
	s.B = b

	s.ready = false
	s.second = 0.0
	s.Med = 0.0
	s.D = 1.0

	s.hx = nil
	s.hy = nil
}
func (s *Autoregr) next() float64 {
	if s.ready {
	   s.ready = false;
	   return s.second * s.D + s.Med
	} else {
	   u := 0.0
	   v := 0.0
	   ss := 0.0
	   for (ss > 1.0 || ss == 0.0) {
		   u = 2.0 * r1.Float64() - 1.0;
		   v = 2.0 * r1.Float64() - 1.0;
		   ss = u * u + v * v;
	   }
	   
	   r := math.Sqrt(-2.0 * math.Log(ss) / ss);
	   s.second = r * u;
	   s.ready = true;
	   return r * v * s.D + s.Med;
	}
}
func (s *Autoregr) GetMax() float64 {
    return 0
}
func (s *Autoregr) GetMin() float64 {
    return 0
}
func (s *Autoregr) Run(in float64) float64 {
	_y := s.next()
	//_y := r1.Float64()
	
	s.hx = append(s.hx, _y)
	Q, P := 0.0, 0.0

	//fmt.Println("A", s.A)
	//fmt.Println("B", s.B)
	//fmt.Println("s.hx", s.hx)
	//fmt.Println("s.hy", s.hy)

	for i := 0; i < len(s.B); i+=1 {
		if int(in)-(i+1) >= 0{
			Q += s.B[i] * s.hx[int(in)-(i+1)]
		} else {
			break
		}
	}
	
	for i := 0; i < len(s.A); i+=1 {
		if int(in)-(i+1) >= 0 {
			P += s.A[i] * s.hy[int(in)-(i+1)]
		} else {
			break
		}
	}
	_y += (Q+P)
	s.hy = append(s.hy, _y) 
	return _y
}

func Gen(start, end float64, n int, sig Signal, s, e float64) (plotter.XYs, []float64) {
	//var s1 TonalEnvelope

	//s1.SetParams(1, 1.5, 0.5, 0.7, 0.7)
	//s1.SetPeriod(2)

	//s1.SetAmplitude(2)
	//s1.SetFrequency(1.5)
	//s1.SetInitialPhase(1.2)
	//s1.SetDelay(delay)
	//s1.SetA(0.5)
	var X []float64

	t := (end-start)/float64(n)

	pts := make(plotter.XYs, int(e-s)*3)
	j := 0
	for i := 0; i < n*3.0; i+=3.0 {
		
		//fmt.Println("id", i, j)

		__x := float64(j)*t
		__y := sig.Run(__x)
		if (j < int((e-s))) {
			pts[i].X = __x
			pts[i].Y = 0
			pts[i+1].X = float64(j)*t
			pts[i+1].Y = __y
		}
		X = append(X, __y)
		if (j < int((e-s))) {
			pts[i+2].X = float64(j)*t
			pts[i+2].Y = 0
		}
		j+=1
	}
	return pts, X
}

func Gen2(start, end float64, n int, sig Signal, s,e float64) plotter.XYs {
	
	t := (end-start)/float64(n)

	pts := make(plotter.XYs, int(e-s)*2-1)
	j := 0
	for i := 0; i < n*2; i+=1 {
		if (j >= int((e-s))) {
			break
		}
		if i == 0 {
			pts[i].X = float64(j)*t
			pts[i].Y = sig.Run(pts[i].X)
			j+=1
			continue
		}
		i+=1
		pts[i].X = float64(j)*t
		pts[i].Y = sig.Run(pts[i].X)

		pts[i-1].X = (pts[i].X- pts[i-2].X)/2 + pts[i-2].X
		pts[i-1].Y = pts[i-2].Y+(pts[i].Y- pts[i-2].Y) /
			(pts[i].X - pts[i-2].X) * (pts[i-1].X - pts[i-2].X)

		j+=1
	}
	//for i:= 0; i < 2*n+1; i +=1 {
	//	fmt.Println(pts[i].X, pts[i].Y)
	//}
	return pts
}

func GenForFFT(start, end float64, n int, data []complex128, s, e float64) (plotter.XYs,plotter.XYs) {
	t := (end-start)/float64(n)

	ptsampl := make(plotter.XYs, int(e-s))
	ptsphas := make(plotter.XYs, int(e-s))

	j := 0
	for i := 0; i < n; i+=1 {
		if (j >= int((e-s))) {
			break
		}
		ptsampl[i].X = float64(j)*t
		ptsampl[i].Y = cmplx.Abs(data[i])
		
		ptsphas[i].X = float64(j)*t
		ptsphas[i].Y = math.Atan(imag(data[i])/real(data[i]))
		j+=1
	}

	return ptsampl, ptsphas
}

func GenForDens(start, end float64, data, fr []float64, s, e float64) (plotter.XYs) {
	//fmt.Println("@@: ", data)
	res := make(plotter.XYs, int(e-s))
	
	t := (end-start)/float64(len(data))
	j := 0
	for i := 0; i < len(data); i+=1 {
		if (j >= int((e-s))) {
			break
		}
		res[i].X = float64(j)*t
		res[i].Y = data[i]
		j+=1
	}

	return res
}

func smooth(input []complex128, n int, window int) []complex128 {
   var z,k1,k2,hw float64
   var tmp complex128
   if math.Mod(float64(window), 2.0) == 0.0 {
	   window+=1
   }
   
   hw = float64(window-1) / 2;
   
   var output []complex128

   output = append(output, input[0]);

   for i := 1; i < n; i+=1 {
       tmp = complex(0, 0)
       if float64(i) < hw {
           k1 = 0
           k2 = 2 * float64(i)
           z = k2 + 1
       } else if (float64(i)+hw) > (float64(n)-1.0) {
           k1 = float64(i)-float64(n)+float64(i)+1.0
           k2 = float64(n)-1.0
           z = k2-k1+1.0
       } else {
           k1 = float64(i)-hw
           k2 = float64(i)+hw
           z = float64(window)
       }

       for j := int(k1); j <= int(k2); j+=1 {
           tmp = tmp + input[j]
       }
       output = append(output, tmp / complex(z, 0))
   }

   return output
}

type Lab4 struct {
	N int
	X []float64
	Choice int
	A, Ampl, Freq, Phase, Med, D, Delay, Period float64
	EnvelopeWidth float64
	EnvelopeFreq float64
	IA, IB float64
	M float64
	
	ArrA []float64
	ArrB []float64

	Start, End float64

	DelayedSingleImpulse_ DelayedSingleImpulse
	DelayedUnitJump_ DelayedUnitJump
	DiscretizedDecreasingExponent_ DiscretizedDecreasingExponent
	SampledSine_ SampledSine
	Meander_ Meander
	Saw_ Saw
	ExponentialEnvelope_ ExponentialEnvelope
	BalancedEnvelope_ BalancedEnvelope
	TonalEnvelope_ TonalEnvelope
	WhiteNoise_ WhiteNoise
	WhiteNoiseNorm_ WhiteNoiseNorm
	Autoregr_ Autoregr
	ParamMiddle, ParamD, ParamSigma, ParamVar, ParamSim, ParamEks float64

	ImStart, ImEnd float64
	Smooth int
}
func (in *Lab4) Init() {
	in.N = 0
	in.Choice = 0

	in.Start = 0
	in.End = 1
	in.N = 10
}

func (in *Lab4) Go() {

	var s Signal 

	switch in.Choice {
	case 1:
		in.DelayedSingleImpulse_.SetParams(in.Delay)
		s = &(in.DelayedSingleImpulse_)
	case 2:
		in.DelayedUnitJump_.SetParams(in.Delay)
		s = &(in.DelayedUnitJump_)
	case 3:
		in.DiscretizedDecreasingExponent_.SetParams(in.A)
		s = &(in.DiscretizedDecreasingExponent_)
	case 4:
		in.SampledSine_.SetParams(in.Ampl, in.Freq, in.Phase)
		s = &(in.SampledSine_)
	case 5:
		in.Meander_.SetParams(in.Period)
		s = &(in.Meander_)
	case 6:
		in.Saw_.SetParams(in.Period)
		s = &(in.Saw_)
	case 7:
		in.ExponentialEnvelope_.SetParams(in.Ampl, 
			in.Freq, in.Phase, in.EnvelopeWidth)
		s = &(in.ExponentialEnvelope_)
	case 8:
		in.BalancedEnvelope_.SetParams(in.Ampl, 
			in.Freq, in.Phase, in.EnvelopeFreq)
		s = &(in.BalancedEnvelope_)
	case 9:
		in.TonalEnvelope_.SetParams(in.Ampl,
			in.Freq, in.Phase, in.EnvelopeFreq, in.M)
		s = &(in.TonalEnvelope_)

	case 10:
		in.WhiteNoise_.SetParams(in.IA, in.IB)
		s = &(in.WhiteNoise_)
	case 11:
		in.WhiteNoiseNorm_.SetParams(in.Med, in.D)
		s = &(in.WhiteNoiseNorm_)
	case 12:
		in.Autoregr_.SetParams(in.ArrA, in.ArrB)
		s = &(in.Autoregr_)
	}

	var lineData plotter.XYs
	lineData, in.X = Gen(in.Start, in.End, in.N, s, in.ImStart, in.ImEnd)
	GenImage("Digital signal impulse", "resources/data1.png", lineData)
	lineDataLine := Gen2(in.Start, in.End, in.N, s, in.ImStart, in.ImEnd)
	GenImage("Digital signal line","resources/data2.png", lineDataLine)

	//Some params
	val := 0.0
	for i:=0; i < len(in.X); i+=1 {
		val += in.X[i]
	}
	in.ParamMiddle = val / float64(len(in.X))
	
	val = 0.0
	for i:=0; i < len(in.X); i+=1 {
		val += (in.X[i] - in.ParamMiddle)*(in.X[i] - in.ParamMiddle)
	}
	in.ParamD = val / float64((len(in.X)-1))

	in.ParamSigma = math.Sqrt(in.ParamD)

	in.ParamVar = in.ParamSigma / in.ParamMiddle

	val = 0.0
	for i:=0; i < len(in.X); i+=1 {
		val += (in.X[i] - in.ParamMiddle)*
			(in.X[i] - in.ParamMiddle)*(in.X[i] - in.ParamMiddle)
	}
	val /= float64(len(in.X))
	in.ParamSim = val / (in.ParamSigma*in.ParamSigma*in.ParamSigma)

	val = 0.0
	for i:=0; i < len(in.X); i+=1 {
		val += (in.X[i] - in.ParamMiddle)*(in.X[i] - in.ParamMiddle)*
			(in.X[i] - in.ParamMiddle)*(in.X[i] - in.ParamMiddle)
	}
	val /= float64(len(in.X))
	val /= (in.ParamSigma*in.ParamSigma*in.ParamSigma*in.ParamSigma)
	in.ParamEks = val - 3.0
	//end params

	var fftarr []complex128
	for i := 0; i < len(in.X); i+=1 {
		fftarr = append(fftarr, complex(in.X[i], 0))
	}
	fftres := fft.FFT(fftarr)
	//fmt.Println("FFT")
	//fmt.Println(fftres)
	ampl, phas := GenForFFT(in.Start, in.End, in.N, fftres, in.ImStart, in.ImEnd)
	GenImage("FFT amplitude","resources/data3.png", ampl)
	GenImage("FFT phase","resources/data4.png", phas)

	fftreswin := smooth(fftres, len(fftres), in.Smooth)

	amplwin, phaswin := GenForFFT(in.Start, in.End, in.N, fftreswin, in.ImStart, in.ImEnd)
	GenImage("FFT amplitude with smooth","resources/data5.png", amplwin)
	GenImage("FFT phase with smooth","resources/data6.png", phaswin)

	//t := (in.End-in.Start)/float64(in.N)
	var po spectral.PwelchOptions
	pxx, freqs := spectral.Pwelch(in.X, (in.End-in.Start)/float64(in.N), &po)
	

	GenImage("Pwelch","resources/data7.png", GenForDens(in.Start, in.End, pxx, freqs, in.ImStart, in.ImEnd))
	//GenImage("resources/data8.png", GenForDens(in.Start, in.End, in.N, freqs))
}

func GenImage(title string, name string, data plotter.XYs) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = title
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	l, err := plotter.NewLine(data)
	if err != nil {
		panic(err)
	}
	l.LineStyle.Width = vg.Points(1)
	l.LineStyle.Color = color.RGBA{B: 255, A: 255}
	p.Add(l)
	p.Legend.Add("line", l)
	if err := p.Save(4*vg.Inch, 4*vg.Inch, name); err != nil {
		panic(err)
	}
}

/*
func main() {

	var s1 TonalEnvelope
	s1.SetParams(1, 1.5, 0.5, 0.7, 0.7)

	lineData := Gen2(20, 5, 10, &s1)

	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	p.Title.Text = "Digital signal"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"

	l, err := plotter.NewLine(lineData)
	if err != nil {
		panic(err)
	}
	l.LineStyle.Width = vg.Points(1)
	//l.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
	l.LineStyle.Color = color.RGBA{B: 255, A: 255}

	p.Add(l)
	p.Legend.Add("line", l)

	// Save the plot to a PNG file.
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "points.png"); err != nil {
		panic(err)
	}

	n:=20
	var fftarr []complex128
	j:=0
	t:=0.5
	for i := 0; i < n; i+=1 {
		//fmt.Println(s1.Run(float64(j)*t/float64(n)))
		fftarr = append(fftarr, complex(s1.Run(float64(j)*t/float64(n)), 0))
		j+=1
	}
	fftres := fft.FFT(fftarr)


	//fftcontrol := fft.IFFT(fftres)
	//fmt.Println("!!!!!!!!!")
	//fmt.Println(fftcontrol)
}*/
