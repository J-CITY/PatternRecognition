package main

import (
	//"fmt"
	"os"
	"sort"
	"log"
	"github.com/google/gxui"
	"github.com/google/gxui/drivers/gl"
	"github.com/google/gxui/math"
	"github.com/google/gxui/samples/flags"
)
import "path/filepath"

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

func (w *FileManager)Go(in string) {
	if (in == ".") {
		for i := len(w.path)-1; i >=0; i-=1 {
			if w.path[i] == '\\' {
				w.path = w.path[:i]
				break
			}
		}
		w.SetPath(w.path)
		return
	}
	isdir := w.IsDir(w.path+"\\"+in)
	if (isdir) {
		w.SetPath(w.path+"\\"+in)
	} else {

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


func GetButton(text string, theme gxui.Theme) gxui.Button {
	btn := theme.CreateButton()
	btn.SetText(text)
	return btn
}
// Number picker uses the gxui.DefaultAdapter for driving a list
func numberPicker(theme gxui.Theme, overlay gxui.BubbleOverlay) gxui.Control {
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

	layoutControlOpen := theme.CreateLinearLayout()
    layoutControlOpen.SetDirection(gxui.LeftToRight)
	labelOpen := theme.CreateLabel()
	labelOpen.SetText("Open: ")
	layoutControlOpen.AddChild(labelOpen)

	btnOpen := GetButton("OK", theme)
	btnOpen.OnClick(func(gxui.MouseEvent) {
		
	})
	layoutControlOpen.AddChild(btnOpen)

	layout.AddChild(layoutControlOpen)

	list.OnItemClicked(func(ev gxui.MouseEvent, item gxui.AdapterItem) {
		//if dropList.Selected() != item {
		//	dropList.Select(item)
		//}
		if ev.Button == gxui.MouseButtonRight {
			fm.Go(item.(string))
			adapter.SetItems(fm.GetList())
			
		}
		
		//selected.SetText(fmt.Sprintf("%s - %d", item, adapter.ItemIndex(item)))
	})

	return layout
}



func appMain(driver gxui.Driver) {
	theme := flags.CreateTheme(driver)

	overlay := theme.CreateBubbleOverlay()

	

	holder := theme.CreatePanelHolder()
	holder.AddPanel(numberPicker(theme, overlay), "Default adapter")

	window := theme.CreateWindow(800, 600, "Lists")
	window.SetScale(flags.DefaultScaleFactor)
	window.AddChild(holder)
	window.AddChild(overlay)
	window.OnClose(driver.Terminate)
	window.SetPadding(math.Spacing{L: 10, T: 10, R: 10, B: 10})
}

func main() {
	gl.StartDriver(appMain)
}