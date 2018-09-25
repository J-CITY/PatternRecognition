package UI

import (
	"github.com/google/gxui"
	"strconv"
	"fmt"
	"os"
	//"log"
	"io/ioutil"
	"sort"
	"log"
	"io"
	"path/filepath"
	"strings"
)

import (
	"image/color"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	//"gonum.org/v1/plot/plotutil"
	"image"
	"image/draw"
	"gonum.org/v1/plot/vg"
	"lab/Lab2No1"
	"lab/Lab2No2"
	"lab/Lab4"
)
//FILE MANAGER
type FileManager struct {
	list []string
	path, RootPath string
}

func (w *FileManager)GetList() []string {
	return w.list
}

func (w *FileManager)SetPath(_path string) {
	w.path = _path

	files, _ := w.ReadDir(w.path)
	w.list = nil
	if w.RootPath != w.path {
		w.list = append(w.list, ".")
	}
	for i := 0; i < len(files); i+=1 {
		w.list = append(w.list, files[i].Name())
	}
}

func (w *FileManager)IsDir(path string) bool {
	finfo, err := os.Stat(path)
	if err != nil {
	    return false
	}
	if finfo.IsDir() {
		return true
	} else {
	    return false
	}
}

func (w *FileManager)Go(in string) bool {
	if (in == ".") {
		for i := len(w.path)-1; i >=0; i-=1 {
			if w.path[i] == '\\' {
				w.path = w.path[:i]
				break
			}
		}
		w.SetPath(w.path)
		return true
	}
	isdir := w.IsDir(w.path+"\\"+in)
	if (isdir) {
		w.SetPath(w.path+"\\"+in)
		return true
	} else {
		return false
	}
}

func (w *FileManager)Init() {
	dir, err := filepath.Abs(filepath.Dir(os.Args[0]))
    if err != nil {
        log.Fatal(err)
	}
	w.RootPath = dir
}

func (w *FileManager)ReadDir(dirname string) ([]os.FileInfo, error) {
    f, err := os.Open(dirname)
    if err != nil {
        return nil, err
    }
    list, err := f.Readdir(-1)
    f.Close()
    if err != nil {
        return nil, err
    }
    sort.Slice(list, func(i, j int) bool { return list[i].Name() < list[j].Name() })
    return list, nil
}

func FileList(theme gxui.Theme, pathFileData *string) gxui.Control {
	var fm FileManager
	fm.Init()
	fm.SetPath(fm.RootPath)

	adapter := gxui.CreateDefaultAdapter()
	adapter.SetItems(fm.GetList())

	layout := theme.CreateLinearLayout()
	layout.SetDirection(gxui.TopToBottom)

	layoutControl := theme.CreateLinearLayout()
    layoutControl.SetDirection(gxui.LeftToRight)
	labelPath := theme.CreateLabel()
	labelPath.SetText("Set root path: ")
	layoutControl.AddChild(labelPath)

	inputPath := theme.CreateTextBox()
	layoutControl.AddChild(inputPath)

	btn := GetButton("OK", theme)
	btn.OnClick(func(gxui.MouseEvent) {
		_path := inputPath.Text()
		if fm.IsDir(_path) {
			fm.RootPath = _path
			fm.SetPath(fm.RootPath)
			adapter.SetItems(fm.GetList())
		}
	})
	layoutControl.AddChild(btn)

	layout.AddChild(layoutControl)


	list := theme.CreateList()
	list.SetAdapter(adapter)
	list.SetOrientation(gxui.Vertical)
	layout.AddChild(list)

	

	list.OnItemClicked(func(ev gxui.MouseEvent, item gxui.AdapterItem) {
		if ev.Button == gxui.MouseButtonLeft {
			b := fm.Go(item.(string))
			if b {
				adapter.SetItems(fm.GetList())
			} else {
				str := item.(string)
				if str[len(str)-3:] == ".dt" {
					*pathFileData = fm.path + "\\" + str
				}
			}
			
		}
	})

	return layout
}

//UI
func GetLable(text string, theme gxui.Theme) gxui.Label {
	label := theme.CreateLabel()
	label.SetText(text)
	return label
}

func GetButton(text string, theme gxui.Theme) gxui.Button {
	btn := theme.CreateButton()
	btn.SetText(text)
	return btn
}
func GetImage(path string, driver gxui.Driver, theme gxui.Theme) gxui.Image {
	f, err := os.Open(path)
	if err != nil {
		fmt.Printf("Failed to open image '%s': %v\n", path, err)
		os.Exit(1)
	}

	source, _, err := image.Decode(f)
	if err != nil {
		fmt.Printf("Failed to read image '%s': %v\n", path, err)
		os.Exit(1)
	}
	img := theme.CreateImage()
	//img.SetBorderPen(gxui.CreatePen(2.0, gxui.White))

	rgba := image.NewRGBA(source.Bounds())
	draw.Draw(rgba, source.Bounds(), source, image.ZP, draw.Src)
	texture := driver.CreateTexture(rgba, 1)
	img.SetTexture(texture)
	return img
}
func GenHistogram(data []int, N int, HistSize int, fname string, xCoeff float64, min float64, max float64) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "Distribution histogram"
	//p.BackgroundColor = color.RGBA{R:25, G: 25, B: 25, A: 255}
	p.X.Label.Text = "X"
	p.Y.Label.Text = "P"

	xys := make(plotter.XYs, HistSize)
	b := false
	start := 0
	for j := range xys {
		if min != 0 || max != 0 {
			if (float64(j) * xCoeff < min) || (float64(j) * xCoeff > max) {
				continue
			} else if !b {
				b = true
				start = j
			} 
		}
		xys[j].X = float64(j) * xCoeff
		xys[j].Y = float64(data[j]) / float64(N)
	}
	h, err := plotter.NewHistogram(xys[start:], HistSize-start)
	if err != nil {
		panic(err)
	}
	h.LineStyle.Color = color.RGBA{B: 255, A: 255}
	//h.Normalize(1)
	p.Add(h)

	if err := p.Save(4*vg.Inch, 4*vg.Inch, "resources/"+fname); err != nil {
		panic(err)
	}
}
func GetDistributionFunctions(data []int, N int, FuncSize int, fname string, xCoeff float64, min float64, max float64) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "Distribution Functions"
	//p.BackgroundColor = color.RGBA{R:25, G: 25, B: 25, A: 255}
	p.X.Label.Text = "X"
	p.Y.Label.Text = "F(X)"
	p.Add(plotter.NewGrid())

	xys := make(plotter.XYs, 2)
	for i := -1; i < FuncSize+1; i+=1 {
		if min != 0 || max != 0 {
			if (float64(i) * xCoeff < min) || (float64(i) * xCoeff > max) {
				continue
			} 
		}
		if i < 0 {
			xys[0].X = float64(i) * xCoeff
			xys[0].Y = 0
			xys[1].X = float64(i+1) * xCoeff
			xys[1].Y = 0
		} else if i == FuncSize+1 {
			xys[0].X = float64(i) * xCoeff
			xys[0].Y = 1
			xys[1].X = float64(i+1) * xCoeff
			xys[1].Y = 1
		} else {
			xys[0].X = float64(i) * xCoeff
			xys[0].Y += float64(data[i]) / float64(N)
			
			xys[1].X = float64(i+1) * xCoeff
			xys[1].Y += float64(data[i]) / float64(N)
		}
		
		l, err := plotter.NewLine(xys)
		if err != nil {
			panic(err)
		}
		l.LineStyle.Width = vg.Points(1)
		//l.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
		l.LineStyle.Color = color.RGBA{B: 255, A: 255}
		p.Add(l)

	}
	

	if err := p.Save(4*vg.Inch, 4*vg.Inch, "resources/"+fname); err != nil {
		panic(err)
	}
}

func AddLab2No1(driver gxui.Driver, theme gxui.Theme) gxui.LinearLayout {
	var logic Lab2No1.Lab2No1
	logic.Init()

	layoutSplit := theme.CreateLinearLayout()
	layoutSplit.SetDirection(gxui.LeftToRight)

	layoutScroll := theme.CreateScrollLayout() 

	layoutMain := theme.CreateLinearLayout()
	layoutInfo := theme.CreateLinearLayout()
	layoutHist := theme.CreateLinearLayout()
	layoutDistFunc := theme.CreateLinearLayout()
	
	layoutMain.SetDirection(gxui.TopToBottom)

	layoutControl := theme.CreateLinearLayout()
    layoutControl.SetDirection(gxui.LeftToRight)
	layoutControl.AddChild(GetLable("Set N: ", theme))
	inputN := theme.CreateTextBox()
	layoutControl.AddChild(inputN)
	btn := GetButton("OK", theme)
	
	updateImgs := func() {
		GenHistogram(logic.Frequencies, logic.N, 30, "histLab2No1.png", 1, logic.HXmin, logic.HXmax)
		layoutHist.AddChildAt(1, GetImage("resources/histLab2No1.png", driver, theme))
		layoutHist.RemoveChildAt(0)
		GetDistributionFunctions(logic.Frequencies, logic.N, 30, "distFuncLab2No1.png", 1, logic.FXmin, logic.FXmax)
		layoutDistFunc.AddChildAt(1, GetImage("resources/distFuncLab2No1.png", driver, theme))
		layoutDistFunc.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.AddChild(GetLable("M(x) = "    + strconv.FormatFloat(logic.MiddleX, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("D(x) = "         + strconv.FormatFloat(logic.D, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("σ(x) = "         + strconv.FormatFloat(logic.Sigma, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("γ = " + strconv.FormatFloat(logic.Asymmetry, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("k = "  + strconv.FormatFloat(logic.Kurtosis, 'f', 6, 64), theme))
	}

	btn.OnClick(func(gxui.MouseEvent) {
		if (inputN.Text() != "") {
			logic.N, _ = strconv.Atoi(inputN.Text())
			fmt.Print("N: ", logic.N)
			logic.GenX()
			
			updateImgs()
		}
	})
	layoutControl.AddChild(btn)
	layoutMain.AddChild(layoutControl)

	layoutInfo.SetHorizontalAlignment(gxui.AlignLeft)
	layoutInfo.AddChild(GetLable("M(x) = ", theme))
	layoutInfo.AddChild(GetLable("D(x) = ", theme))
	layoutInfo.AddChild(GetLable("σ(x) = ", theme))
	layoutInfo.AddChild(GetLable("γ = ", theme))
	layoutInfo.AddChild(GetLable("k = ", theme))
	layoutMain.AddChild(layoutInfo)

	//Histo range
	layoutHistControl := theme.CreateLinearLayout()
    layoutHistControl.SetDirection(gxui.LeftToRight)
	layoutHistControl.AddChild(GetLable("Set range: ", theme))
	inputHX0 := theme.CreateTextBox()
	layoutHistControl.AddChild(inputHX0)
	inputHXN := theme.CreateTextBox()
	layoutHistControl.AddChild(inputHXN)
	btnHR := GetButton("OK", theme)
	btnHR.OnClick(func(gxui.MouseEvent) {
		if (inputHX0.Text() != "") || (inputHXN.Text() != "") {
			logic.HXmin = 0
			logic.HXmax = 0
		}
		if (inputHX0.Text() != "") {
			logic.HXmin, _ = strconv.ParseFloat(inputHX0.Text(), 64)
		}
		if (inputHXN.Text() != "") {
			logic.HXmax, _ = strconv.ParseFloat(inputHXN.Text(), 64)
		}
	})
	layoutHistControl.AddChild(btnHR)
	layoutMain.AddChild(layoutHistControl)

	layoutHist.SetHorizontalAlignment(gxui.AlignCenter)
	layoutHist.AddChild(GetImage("resources/default.png", driver, theme))
	layoutMain.AddChild(layoutHist)

	//Func range
	layoutFuncControl := theme.CreateLinearLayout()
    layoutFuncControl.SetDirection(gxui.LeftToRight)
	layoutFuncControl.AddChild(GetLable("Set range: ", theme))
	inputFX0 := theme.CreateTextBox()
	layoutFuncControl.AddChild(inputFX0)
	inputFXN := theme.CreateTextBox()
	layoutFuncControl.AddChild(inputFXN)
	btnFR := GetButton("OK", theme)
	btnFR.OnClick(func(gxui.MouseEvent) {
		if (inputFX0.Text() == "") || (inputFXN.Text() == "") {
			logic.FXmin = 0
			logic.FXmax = 0
		}

		if (inputFX0.Text() != "") {
			logic.FXmin, _ = strconv.ParseFloat(inputFX0.Text(), 64)
		}
		if (inputFXN.Text() != "") {
			logic.FXmax, _ = strconv.ParseFloat(inputFXN.Text(), 64)
		}
	})
	layoutFuncControl.AddChild(btnFR)
	layoutMain.AddChild(layoutFuncControl)

	layoutDistFunc.SetHorizontalAlignment(gxui.AlignCenter)
	layoutDistFunc.AddChild(GetImage("resources/default.png", driver, theme))
	layoutMain.AddChild(layoutDistFunc)

	pathFileData := ""

	layoutList := theme.CreateLinearLayout()
    layoutList.SetDirection(gxui.TopToBottom)
	layoutList.AddChild(FileList(theme, &pathFileData))

	layoutControlOpen := theme.CreateLinearLayout()
    layoutControlOpen.SetDirection(gxui.LeftToRight)
	labelOpen := theme.CreateLabel()
	labelOpen.SetText("Open: ")
	layoutControlOpen.AddChild(labelOpen)

	btnOpen := GetButton("OK", theme)
	btnOpen.OnClick(func(gxui.MouseEvent) {
		fmt.Println(pathFileData)
		X, err := ioutil.ReadFile(pathFileData)
		if err != nil {
			panic(err)
		}

		//fmt.Print("N: ", logic.N)
		
		logic.SetX(X)
		updateImgs()
	})
	layoutControlOpen.AddChild(btnOpen)
	layoutList.AddChild(layoutControlOpen)
	//
	layout := theme.CreateLinearLayout()
    layout.SetDirection(gxui.LeftToRight)
	layout.AddChild(layoutList)
	//
	//layoutSplit.AddChild(layoutList)
	layoutSplit.AddChild(layoutMain)
	layoutScroll.SetChild(layoutSplit)
	//
	layout.AddChild(layoutScroll)
	//
	return layout
}

func AddLab2No2(driver gxui.Driver, theme gxui.Theme) gxui.LinearLayout {
	var logic Lab2No2.Lab2No2
	logic.Init()

	layoutSplit := theme.CreateLinearLayout()
	layoutSplit.SetDirection(gxui.LeftToRight)

	layoutScroll := theme.CreateScrollLayout() 

	layoutMain := theme.CreateLinearLayout()
	layoutInfo := theme.CreateLinearLayout()
	layoutHist := theme.CreateLinearLayout()
	layoutDistFunc := theme.CreateLinearLayout()
	
	layoutMain.SetDirection(gxui.TopToBottom)

	layoutControl := theme.CreateLinearLayout()
    layoutControl.SetDirection(gxui.LeftToRight)
	layoutControl.AddChild(GetLable("Set N (Размер выборки): ", theme))
	inputN := theme.CreateTextBox()
	layoutControl.AddChild(inputN)

	layoutControl.AddChild(GetLable("Set histogrem size: ", theme))
	inputHN := theme.CreateTextBox()
	layoutControl.AddChild(inputHN)

	btn := GetButton("OK", theme)
	
	updateImgs := func() {
		GenHistogram(logic.Frequencies, logic.N, logic.FN-1, "histLab2No2.png", (logic.GetXN()-logic.GetX0())/float64(logic.FN), logic.HXmin, logic.HXmax)
		layoutHist.AddChildAt(1, GetImage("resources/histLab2No2.png", driver, theme))
		layoutHist.RemoveChildAt(0)
		GetDistributionFunctions(logic.Frequencies, logic.N, logic.FN-1, "distFuncLab2No2.png", (logic.GetXN()-logic.GetX0())/float64(logic.FN), logic.FXmin, logic.FXmax)
		layoutDistFunc.AddChildAt(1, GetImage("resources/distFuncLab2No2.png", driver, theme))
		layoutDistFunc.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.RemoveChildAt(0)
		layoutInfo.AddChild(GetLable("M(X) = "    + strconv.FormatFloat(logic.MiddleX, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("D(X) = "         + strconv.FormatFloat(logic.D, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("σ(X) = "         + strconv.FormatFloat(logic.Sigma, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("γ = " + strconv.FormatFloat(logic.Asymmetry, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("k = "  + strconv.FormatFloat(logic.Kurtosis, 'f', 6, 64), theme))
	}

	btn.OnClick(func(gxui.MouseEvent) {
		if (inputN.Text() != "") {
			logic.N, _ = strconv.Atoi(inputN.Text())
			fmt.Print("N: ", logic.N)
			if inputHN.Text() != "" {
				res, _ := strconv.Atoi(inputHN.Text())
				logic.GenX(res)
			} else {
				logic.GenX(0)
			}
			updateImgs()
		}
	})
	layoutControl.AddChild(btn)
	layoutMain.AddChild(layoutControl)

	layoutControlRange := theme.CreateLinearLayout()
    layoutControlRange.SetDirection(gxui.LeftToRight)
	layoutControlRange.AddChild(GetLable("Set range: ", theme))
	layoutControlRange.AddChild(GetLable("x0: ", theme))
	inputX0 := theme.CreateTextBox()
	layoutControlRange.AddChild(inputX0)
	layoutControlRange.AddChild(GetLable("xN: ", theme))
	inputXN := theme.CreateTextBox()
	layoutControlRange.AddChild(inputXN)
	layoutControlRange.AddChild(GetLable("C: ", theme))
	inputC := theme.CreateTextBox()
	layoutControlRange.AddChild(inputC)
	btnRange := GetButton("OK", theme)
	btnRange.OnClick(func(gxui.MouseEvent) {
		if (inputX0.Text() != "") {
			res, _ := strconv.ParseFloat(inputX0.Text(), 64)
			logic.SetX0(res)
		}
		if (inputXN.Text() != "") {
			res, _ := strconv.ParseFloat(inputXN.Text(), 64)
			logic.SetXN(res)
		}
		if (inputC.Text() != "") {
			res, _ := strconv.ParseFloat(inputC.Text(), 64)
			logic.SetC(res)
		}
	})
	layoutControlRange.AddChild(btnRange)
	layoutMain.AddChild(layoutControlRange)


	layoutInfo.SetHorizontalAlignment(gxui.AlignLeft)
	layoutInfo.AddChild(GetLable("M(X) = ", theme))
	layoutInfo.AddChild(GetLable("D(X) = ", theme))
	layoutInfo.AddChild(GetLable("σ(X) = ", theme))
	layoutInfo.AddChild(GetLable("γ =    ", theme))
	layoutInfo.AddChild(GetLable("k =    ", theme))
	layoutMain.AddChild(layoutInfo)

	//Histo range
	layoutHistControl := theme.CreateLinearLayout()
    layoutHistControl.SetDirection(gxui.LeftToRight)
	layoutHistControl.AddChild(GetLable("Set range: ", theme))
	inputHX0 := theme.CreateTextBox()
	layoutHistControl.AddChild(inputHX0)
	inputHXN := theme.CreateTextBox()
	layoutHistControl.AddChild(inputHXN)
	btnHR := GetButton("OK", theme)
	btnHR.OnClick(func(gxui.MouseEvent) {
		if (inputHX0.Text() == "") || (inputHXN.Text() == "") {
			logic.HXmin = 0
			logic.HXmax = 0
		}
		if (inputHX0.Text() != "") {
			logic.HXmin, _ = strconv.ParseFloat(inputHX0.Text(), 64)
		}
		if (inputHXN.Text() != "") {
			logic.HXmax, _ = strconv.ParseFloat(inputHXN.Text(), 64)
		}
	})
	layoutHistControl.AddChild(btnHR)
	layoutMain.AddChild(layoutHistControl)

	layoutHist.SetHorizontalAlignment(gxui.AlignCenter)
	layoutHist.AddChild(GetImage("resources/default.png", driver, theme))
	layoutMain.AddChild(layoutHist)

	//Func range
	layoutFuncControl := theme.CreateLinearLayout()
    layoutFuncControl.SetDirection(gxui.LeftToRight)
	layoutFuncControl.AddChild(GetLable("Set range: ", theme))
	inputFX0 := theme.CreateTextBox()
	layoutFuncControl.AddChild(inputFX0)
	inputFXN := theme.CreateTextBox()
	layoutFuncControl.AddChild(inputFXN)
	btnFR := GetButton("OK", theme)
	btnFR.OnClick(func(gxui.MouseEvent) {
		if (inputFX0.Text() != "") || (inputFXN.Text() != "") {
			logic.FXmin = 0
			logic.FXmax = 0
		}

		if (inputFX0.Text() != "") {
			logic.FXmin, _ = strconv.ParseFloat(inputFX0.Text(), 64)
		}
		if (inputFXN.Text() != "") {
			logic.FXmax, _ = strconv.ParseFloat(inputFXN.Text(), 64)
		}
	})
	layoutFuncControl.AddChild(btnFR)
	layoutMain.AddChild(layoutFuncControl)

	layoutDistFunc.SetHorizontalAlignment(gxui.AlignCenter)
	layoutDistFunc.AddChild(GetImage("resources/default.png", driver, theme))
	layoutMain.AddChild(layoutDistFunc)

	pathFileData := ""

	layoutList := theme.CreateLinearLayout()
    layoutList.SetDirection(gxui.TopToBottom)
	layoutList.AddChild(FileList(theme, &pathFileData))

	layoutControlOpen := theme.CreateLinearLayout()
    layoutControlOpen.SetDirection(gxui.LeftToRight)
	labelOpen := theme.CreateLabel()
	labelOpen.SetText("Open: ")
	layoutControlOpen.AddChild(labelOpen)

	btnOpen := GetButton("OK", theme)
	btnOpen.OnClick(func(gxui.MouseEvent) {
		fmt.Println(pathFileData)
		file, err := os.Open(pathFileData)
    	if err != nil {
    	    os.Exit(1)
    	}
		defer file.Close()
		
		var perline float64
        var nums []float64

        for {
           _, err := fmt.Fscan(file, &perline) // give a patter to scan
           if err != nil {
              if err == io.EOF {
                      break // stop reading the file
              }
              fmt.Println(err)
              os.Exit(1)
           }
           nums = append(nums, perline)
		}
		
		//_N := int(nums[0])
		_HN := int(nums[1])
		_x0 := nums[2]
		_xN := nums[3]
		_c := nums[4]
		X := nums[5:]
		logic.SetX0(_x0)
		logic.SetXN(_xN)
		logic.SetC(_c)
		logic.SetX(X, _HN)
		updateImgs()
	})
	layoutControlOpen.AddChild(btnOpen)
	layoutList.AddChild(layoutControlOpen)
	layout := theme.CreateLinearLayout()
    layout.SetDirection(gxui.LeftToRight)
	layout.AddChild(layoutList)
	layoutSplit.AddChild(layoutMain)
	layoutScroll.SetChild(layoutSplit)
	layout.AddChild(layoutScroll)

	return layout
}

const s1 string = "01.	задержанный единичный импульс"
const s2 string = "02.	задержанный единичный скачок"
const s3 string = "03.	дискретизированная убывающая экспонента"
const s4 string = "04. дискретизированная синусоида "
const s5 string = "05. меандр"
const s6 string = "06. пила"
const s7 string = "07. сигнал с экспоненциальной огибающей"
const s8 string = "08. cигнал с балансной огибающей"
const s9 string = "09. cигнал с тональной огибающей"
const s10 string = "10. сигнал белого шума"
const s11 string = "11. сигнал белого шума с нормальным распределением"
const s12 string = "12. случайный сигнал авторегрессии-скользящего среднего"

func AddDropList(driver gxui.Driver, theme gxui.Theme, overlay gxui.BubbleOverlay) gxui.DropDownList {
	items := []string{
		s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,
	}

	adapter := gxui.CreateDefaultAdapter()
	adapter.SetItems(items)

	dropList := theme.CreateDropDownList()
	dropList.SetAdapter(adapter)
	dropList.SetBubbleOverlay(overlay)
	return dropList
}

func AddLab4(driver gxui.Driver, theme gxui.Theme, overlay gxui.BubbleOverlay) gxui.LinearLayout {
	var logic Lab4.Lab4
	logic.Init()
	
	choice := 0
	//BUTTNOS ANDLABLES FOR CONTROL LAYOUT
	lableDelay := GetLable("Delay", theme)
	lableA := GetLable("A", theme)
	lableAmpl := GetLable("Amplitude", theme)
	lableFreq := GetLable("Frequency", theme)
	lablePhase := GetLable("Phase", theme)
	lablePeriod := GetLable("Period", theme)
	lableWidth := GetLable("Envelope width", theme)
	lableCarrF := GetLable("Carrier frequency", theme)
	lableEnvF := GetLable("Envelope frequency", theme)
	lableM := GetLable("Modulation depth", theme)
	lableIA := GetLable("A", theme)
	lableIB := GetLable("B", theme)
	lableD := GetLable("D", theme)
	input1 := theme.CreateTextBox()
	input2 := theme.CreateTextBox()
	input3 := theme.CreateTextBox()
	input4 := theme.CreateTextBox()
	input5 := theme.CreateTextBox()

	inputS := theme.CreateTextBox()
	inputE := theme.CreateTextBox()
	inputSmooth := theme.CreateTextBox()

	layoutData1 := theme.CreateLinearLayout()
	layoutData1.SetDirection(gxui.TopToBottom)
	layoutData1.AddChild(GetImage("resources/default.png", driver, theme))

	layoutData2 := theme.CreateLinearLayout()
	layoutData2.SetDirection(gxui.TopToBottom)
	layoutData2.AddChild(GetImage("resources/default.png", driver, theme))

	layoutData3 := theme.CreateLinearLayout()
	layoutData3.SetDirection(gxui.TopToBottom)
	layoutData3.AddChild(GetImage("resources/default.png", driver, theme))
	layoutData3.AddChild(GetImage("resources/default.png", driver, theme))

	layoutData4 := theme.CreateLinearLayout()
	layoutData4.SetDirection(gxui.TopToBottom)
	layoutData4.AddChild(GetImage("resources/default.png", driver, theme))
	layoutData4.AddChild(GetImage("resources/default.png", driver, theme))

	layoutInterval := theme.CreateLinearLayout()
	layoutInterval.SetDirection(gxui.LeftToRight)
	lableStart := GetLable("End", theme)
	lableEnd := GetLable("Start", theme)
	lableN := GetLable("N", theme)
	inputStart := theme.CreateTextBox()
	inputEnd := theme.CreateTextBox()
	inputN := theme.CreateTextBox()
	layoutInterval.AddChild(lableStart)
	layoutInterval.AddChild(inputStart)
	layoutInterval.AddChild(lableEnd)
	layoutInterval.AddChild(inputEnd)
	layoutInterval.AddChild(lableN)
	layoutInterval.AddChild(inputN)

	layoutInfo := theme.CreateLinearLayout()
	layoutInfo.SetDirection(gxui.TopToBottom)
	layoutInfo.AddChild(GetLable("Среднее:", theme))
	layoutInfo.AddChild(GetLable("Дисперсия:", theme))
	layoutInfo.AddChild(GetLable("Среднеквадратичное отклонение:", theme))
	layoutInfo.AddChild(GetLable("Коэффициент вариации:", theme))
	layoutInfo.AddChild(GetLable("Коэффициент асимметрии:", theme))
	layoutInfo.AddChild(GetLable("Коэффициент эксцесса:", theme))

	btnOk := GetButton("OK", theme)
	btnOk.OnClick(func(gxui.MouseEvent) {
		logic.Choice = choice
		r, _ := strconv.ParseFloat(inputStart.Text(), 64)
		logic.Start = r
		r, _ = strconv.ParseFloat(inputEnd.Text(), 64)
		logic.End = r 
		r, _ = strconv.ParseFloat(inputN.Text(), 64)
		logic.N = int(r) 

		switch choice {
		case 1:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Delay = res
		case 2:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Delay = res
		case 3:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.A = res
		case 4:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Ampl = res
			res, _ = strconv.ParseFloat(input2.Text(), 64)
			logic.Freq = res
			res, _ = strconv.ParseFloat(input3.Text(), 64)
			logic.Phase = res
		case 5:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Period = res
		case 6:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Period = res
		case 7:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Ampl = res
			res, _ = strconv.ParseFloat(input2.Text(), 64)
			logic.EnvelopeWidth = res
			res, _ = strconv.ParseFloat(input3.Text(), 64)
			logic.Freq = res
			res, _ = strconv.ParseFloat(input4.Text(), 64)
			logic.Phase = res
		case 8:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Ampl = res
			res, _ = strconv.ParseFloat(input2.Text(), 64)
			logic.EnvelopeFreq = res
			res, _ = strconv.ParseFloat(input3.Text(), 64)
			logic.Freq = res
			res, _ = strconv.ParseFloat(input4.Text(), 64)
			logic.Phase = res
		case 9:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Ampl = res
			res, _ = strconv.ParseFloat(input2.Text(), 64)
			logic.EnvelopeFreq = res
			res, _ = strconv.ParseFloat(input3.Text(), 64)
			logic.Freq = res
			res, _ = strconv.ParseFloat(input4.Text(), 64)
			logic.Phase = res
			res, _ = strconv.ParseFloat(input5.Text(), 64)
			logic.M = res
		case 10:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.IA = res
			res, _ = strconv.ParseFloat(input2.Text(), 64)
			logic.IB = res
		case 11:
			res, _ := strconv.ParseFloat(input1.Text(), 64)
			logic.Med = res
			res, _ = strconv.ParseFloat(input2.Text(), 64)
			logic.D = res
		case 12:
			var arrA []float64
			var arrB []float64
			_A := strings.Split(input1.Text(), " ")
			_B := strings.Split(input2.Text(), " ")

			for i := 0; i < len(_A); i+=1 {
				res, _ := strconv.ParseFloat(_A[i], 64)
				arrA = append(arrA, res)
			}

			for i := 0; i < len(_B); i+=1 {
				res, _ := strconv.ParseFloat(_B[i], 64)
				arrB = append(arrB, res)
			}

			logic.ArrA = arrA
			logic.ArrB = arrB

			//res, _ := strconv.ParseFloat(input3.Text(), 64)
			//logic.D = res
		}

		sm, _ := strconv.ParseFloat(inputSmooth.Text(), 64)
		logic.Smooth = int(sm)
		logic.ImStart, _ = strconv.ParseFloat(inputS.Text(), 64)
		logic.ImEnd, _ = strconv.ParseFloat(inputE.Text(), 64)
		logic.Go()

		layoutData1.RemoveChildAt(0)
		layoutData1.AddChild(GetImage("resources/data1.png", driver, theme))
		layoutData2.RemoveChildAt(0)
		layoutData2.AddChild(GetImage("resources/data2.png", driver, theme))
		layoutData3.RemoveChildAt(0)
		layoutData3.RemoveChildAt(0)
		layoutData3.AddChild(GetImage("resources/data3.png", driver, theme))
		layoutData3.AddChild(GetImage("resources/data5.png", driver, theme))
		layoutData4.RemoveChildAt(0)
		layoutData4.RemoveChildAt(0)
		layoutData4.AddChild(GetImage("resources/data4.png", driver, theme))
		layoutData4.AddChild(GetImage("resources/data6.png", driver, theme))


		for i:=0; i < 6; i+=1 {
			layoutInfo.RemoveChildAt(0)
		}
		layoutInfo.SetDirection(gxui.TopToBottom)
		layoutInfo.AddChild(GetLable("Среднее:"+strconv.FormatFloat(logic.ParamMiddle, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("Дисперсия:"+strconv.FormatFloat(logic.ParamD, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("Среднеквадратичное отклонение:"+strconv.FormatFloat(logic.ParamSigma, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("Коэффициент вариации:"+strconv.FormatFloat(logic.ParamVar, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("Коэффициент асимметрии:"+strconv.FormatFloat(logic.ParamSim, 'f', 6, 64), theme))
		layoutInfo.AddChild(GetLable("Коэффициент эксцесса:"+strconv.FormatFloat(logic.ParamEks, 'f', 6, 64), theme))
	})


	layoutMain := theme.CreateLinearLayout()
	layoutMain.SetDirection(gxui.TopToBottom)

	childSize := 0
	layoutControl := theme.CreateLinearLayout()
	layoutControl.SetDirection(gxui.TopToBottom)

	dropList := AddDropList(driver, theme, overlay)
	dropList.OnSelectionChanged(func(item gxui.AdapterItem) {
		for i := 0; i < childSize; i+=1 {
			layoutControl.RemoveChildAt(0)
		}
		switch item {
		case s1:
			choice = 1
			childSize = 1
			layoutControl.AddChild(lableDelay)
			layoutControl.AddChild(input1)
		case s2:
			choice = 2
			childSize = 1
			layoutControl.AddChild(lableDelay)
			layoutControl.AddChild(input1)
		case s3:
			choice = 3
			childSize = 1
			layoutControl.AddChild(lableA)
			layoutControl.AddChild(input1)
		case s4:
			choice = 4
			childSize = 3
			layoutControl.AddChild(lableAmpl)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableFreq)
			layoutControl.AddChild(input2)
			layoutControl.AddChild(lablePhase)
			layoutControl.AddChild(input3)
		case s5:
			choice = 5
			childSize = 1
			layoutControl.AddChild(lablePeriod)
			layoutControl.AddChild(input1)
		case s6:
			choice = 6
			childSize = 1
			layoutControl.AddChild(lablePeriod)
			layoutControl.AddChild(input1)
		case s7:
			choice = 7
			childSize = 4
			layoutControl.AddChild(lableAmpl)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableWidth)
			layoutControl.AddChild(input2)
			layoutControl.AddChild(lableEnvF)
			layoutControl.AddChild(input3)
			layoutControl.AddChild(lablePhase)
			layoutControl.AddChild(input4)
		case s8:
			choice = 8
			childSize = 4
			layoutControl.AddChild(lableAmpl)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableEnvF)
			layoutControl.AddChild(input2)
			layoutControl.AddChild(lableCarrF)
			layoutControl.AddChild(input3)
			layoutControl.AddChild(lablePhase)
			layoutControl.AddChild(input4)
		case s9:
			choice = 9
			childSize = 5
			layoutControl.AddChild(lableAmpl)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableEnvF)
			layoutControl.AddChild(input2)
			layoutControl.AddChild(lableCarrF)
			layoutControl.AddChild(input3)
			layoutControl.AddChild(lablePhase)
			layoutControl.AddChild(input4)
			layoutControl.AddChild(lableM)
			layoutControl.AddChild(input5)
		case s10:
			choice = 10
			childSize = 2
			layoutControl.AddChild(lableIA)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableIB)
			layoutControl.AddChild(input2)
		case s11:
			choice = 11
			childSize = 2
			layoutControl.AddChild(lableA)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableD)
			layoutControl.AddChild(input2)
		case s12:
			choice = 12
			childSize = 2//4
			layoutControl.AddChild(lableIA)
			layoutControl.AddChild(input1)
			layoutControl.AddChild(lableIB)
			layoutControl.AddChild(input2)
			//layoutControl.AddChild(lableD)
			//layoutControl.AddChild(input3)
		}
		childSize *= 2
		fmt.Println(item)
	})
	layoutMain.AddChild(dropList)
	layoutMain.AddChild(layoutControl)
	layoutMain.AddChild(layoutInterval)
	layoutMain.AddChild(btnOk)

	layoutMain.AddChild(inputSmooth)
	layoutMain.AddChild(inputS)
	layoutMain.AddChild(inputE)


	holder := theme.CreatePanelHolder()
	holder.AddPanel(layoutData1, "DATA1")
	holder.AddPanel(layoutData2, "DATA2")
	holder.AddPanel(layoutInfo,  "INFO")
	holder.AddPanel(layoutData3, "DATA3")
	holder.AddPanel(layoutData4, "DATA4")

	layoutMain.AddChild(holder)

	layoutScroll := theme.CreateScrollLayout() 

	layout := theme.CreateLinearLayout()
    layout.SetDirection(gxui.TopToBottom)

	layoutScroll.SetChild(layoutMain)
	layout.AddChild(layoutScroll)

	return layout
}
