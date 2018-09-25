package main

import (
	"image/color"
	"math/rand"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	//"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"math"
	//"gonum.org/v1/plot/vg/draw"
)


func main() {
	// Get some random points
	rand.Seed(int64(0))
	n := 20
	lineData := randomPoints(n)
	//linePointsData := randomPoints(n)

	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "а(x)"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	p.Add(plotter.NewGrid())


	// Make a line plotter and set its style.
	l, err := plotter.NewLine(lineData)
	if err != nil {
		panic(err)
	}
	l.LineStyle.Width = vg.Points(1)
	l.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
	l.LineStyle.Color = color.RGBA{B: 255, A: 255}

	p.Add(l)
	//p.Legend.Add("line", l)

	// Save the plot to a PNG file.
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "points.png"); err != nil {
		panic(err)
	}
}

// randomPoints returns some random x, y points.
func randomPoints(n int) plotter.XYs {
	step := float64(15) / float64(n)
	
	pts := make(plotter.XYs, n)
	for i := range pts {
		pts[i].X = float64(i)*step
		pts[i].Y = 0.5*pts[i].X*pts[i].X*math.Exp(-pts[i].X)
		//pts[i].Y = -0.5*math.Exp(-pts[i].X)*(pts[i].X*pts[i].X+2*pts[i].X+2)+1
	}
	return pts
}