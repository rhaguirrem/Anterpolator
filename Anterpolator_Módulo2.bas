Attribute VB_Name = "Módulo2"
Sub SortExplorers()
Attribute SortExplorers.VB_ProcData.VB_Invoke_Func = " \n14"
'
' SortExplorers Macro
'

'
    Range("AM4:AP4").Select
    Range(Selection, Selection.End(xlDown)).Select
    ActiveWorkbook.Worksheets("Mapa").Sort.SortFields.Clear
    ActiveWorkbook.Worksheets("Mapa").Sort.SortFields.Add2 Key:=Range("AP4:AP103" _
        ), SortOn:=xlSortOnValues, Order:=xlAscending, DataOption:=xlSortNormal
    With ActiveWorkbook.Worksheets("Mapa").Sort
        .SetRange Range("AM3:AP103")
        .Header = xlYes
        .MatchCase = False
        .Orientation = xlTopToBottom
        .SortMethod = xlPinYin
        .Apply
    End With
End Sub
