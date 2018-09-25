package main

import (
	"github.com/google/gxui"
	"github.com/google/gxui/drivers/gl"
	"flag"
	"github.com/google/gxui/themes/dark"
	"github.com/google/gxui/themes/light"
)

import (
	"lab/UI"
)
	
import "time"
import "math/rand"

var DefaultScaleFactor float32
var FlagTheme string

func CreateTheme(driver gxui.Driver) gxui.Theme {
	if FlagTheme == "light" {
		return light.CreateTheme(driver)
	}
	return dark.CreateTheme(driver)
}

func init() {
	defaultScaleFactor := flag.Float64("scaling", 1.0, "Adjusts the scaling of UI rendering")
	flag.Parse()
	DefaultScaleFactor = float32(*defaultScaleFactor)
	FlagTheme = "dark"
}

func appMain(driver gxui.Driver) {
	theme := CreateTheme(driver)
	overlay := theme.CreateBubbleOverlay()

	holder := theme.CreatePanelHolder()
	
	//holder.AddPanel(UI.AddLab2No1(driver, theme), "Lab2 No1")
	//holder.AddPanel(UI.AddLab2No2(driver, theme), "Lab2 No2")
	holder.AddPanel(UI.AddLab4(driver, theme, overlay), "Lab4")

	window := theme.CreateWindow(800, 600, "Pattern recognition")
	window.SetScale(DefaultScaleFactor)
	window.AddChild(holder)
	window.AddChild(overlay)
	window.OnClose(driver.Terminate)
}

func main() {
	rand.Seed(time.Now().UTC().UnixNano())
	gl.StartDriver(appMain)
}