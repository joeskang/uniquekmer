import PySimpleGUI as sg
import sys
from tqdm import trange


def example1():

    layout = [

        [sg.Text("Enter name, address, phone")],
        [sg.Text('Name', size=(15,1)), sg.InputText('name')],
        [sg.Text('Address', size=(15, 1)), sg.InputText('address')],
        [sg.Text('Phone', size=(15, 1)), sg.InputText('phone')],
        [sg.Submit(), sg.Cancel()]
    ]

    window = sg.Window('Simple data entry window').Layout(layout)
    button, values = window.Read()

    print(button, values[0], values[1], values[2])


def example2():

    if len(sys.argv) == 1:
        event, (fname) = sg.Window("File Opener").Layout([[sg.Text('Documents to open')],
                                                        [sg.Text('Genome File', size=(12, 1)), sg.InputText('~/Documents/uniquekmer/22.fa'), sg.FileBrowse()],
                                                        [sg.Text('Suffix Array', size=(12, 1)), sg.InputText('~/Documents/uniquekmer/22.sa'), sg.FileBrowse()],
                                                        [sg.CloseButton('Open'), sg.CloseButton('Cancel')]]).Read()

    else:
        fname = sys.argv[1]

    if not fname:
        sg.Popup("Cancel", "No filename supplied")
        raise SystemExit("Cancelling: no filename supplied")

    for i in trange(1000):
        sg.OneLineProgressMeter(title="one line",current_value=i+1,max_value=1000,orientation='h', key="progbar")

    print(event, fname)


if __name__ == "__main__":
    example2()
