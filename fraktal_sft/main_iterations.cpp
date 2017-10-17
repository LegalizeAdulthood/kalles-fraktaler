#include "main.h"
#include "main_iterations.h"
#include "fraktal_sft.h"
#include "resource.h"
#include "../formula/formula.h"
#include "../common/tooltip.h"

static int g_nPrevCalc = -1;

static void UpdateIterations(HWND hWnd)
{
	int i = g_SFT.GetIterations();
	SetDlgItemInt(hWnd,IDC_EDIT1,i,FALSE);
}

static void RefreshIterations(HWND hWnd)
{
	char sz[256];
	GetDlgItemText(hWnd,IDC_EDIT1,sz,sizeof(sz));
	int i = atoi(sz);
  g_SFT.SetIterations(i);
}

static void UpdateMinIterations(HWND hWnd, int i)
{
  SetDlgItemInt(hWnd,IDC_EDIT2,i,FALSE);
}

static void UpdateMaxIterations(HWND hWnd, int i)
{
  SetDlgItemInt(hWnd,IDC_EDIT5,i,FALSE);
}

static void UpdateApproxIterations(HWND hWnd)
{

	int i = g_SFT.GetMaxApproximation();
	SetDlgItemInt(hWnd,IDC_EDIT7,i,FALSE);
}

static void UpdateSmoothMethod(HWND hWnd)
{
	int i = g_SFT.GetSmoothMethod();
	SendDlgItemMessage(hWnd,IDC_COMBO2,CB_SETCURSEL,i,0);
}

static void RefreshSmoothMethod(HWND hWnd)
{
	int i = SendDlgItemMessage(hWnd,IDC_COMBO2,CB_GETCURSEL,0,0);
	g_SFT.SetSmoothMethod(i);
}

static void UpdatePower(HWND hWnd)
{
	int i = g_SFT.GetPower();
	SetDlgItemInt(hWnd,IDC_COMBO3,i,FALSE);
}

static int RefreshPower(HWND hWnd, bool refresh = true)
{
	char szPower[256];
	GetDlgItemText(hWnd,IDC_COMBO3,szPower,sizeof(szPower));
	int i = atoi(szPower);
	if (refresh) g_SFT.SetPower(i);
	return i;
}

static void UpdateFractalType(HWND hWnd, int i = -1, int p = -1)
{
	if (i < 0) i = g_SFT.GetFractalType();
	if (p < 0) p = g_SFT.GetPower();
	SendDlgItemMessage(hWnd,IDC_COMBO5,CB_SETCURSEL,i,0);
	int nType = SendDlgItemMessage(hWnd,IDC_COMBO5,CB_GETCURSEL,0,0);
	SendDlgItemMessage(hWnd,IDC_COMBO3,CB_RESETCONTENT,0,0);
	SendDlgItemMessage(hWnd,IDC_COMBO3,CB_ADDSTRING,0,(LPARAM)"2");
	if(nType<=4 || nType>=10)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_ADDSTRING,0,(LPARAM)"3");
	if(nType==0){
		int p;
		for(p=4;p<=10;p++){
			char szNum[4];
			wsprintf(szNum,"%d",p);
			SendDlgItemMessage(hWnd,IDC_COMBO3,CB_ADDSTRING,0,(LPARAM)szNum);
		}
	}
	else{
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_ADDSTRING,0,(LPARAM)"4");
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_ADDSTRING,0,(LPARAM)"5");
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_ADDSTRING,0,(LPARAM)"6");
	}
	SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,p-2,0);
	EnableWindow(GetDlgItem(hWnd,IDC_COMBO3),nType<=4);
	if(nType>=10)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,1,0);
	if(nType>=15)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,2,0);
	if(nType>=27)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,3,0);
	if(nType>=33)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,2,0);
	if(nType==40)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,4,0);
	if(nType==41)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,0,0);
	if(nType>=42)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,1,0);
	if(nType==45 || nType==46)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,2,0);
	if(nType==47 || nType==48)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,3,0);
	if(nType==49 || nType==50)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,4,0);
	if(nType>=51)
		SendDlgItemMessage(hWnd,IDC_COMBO3,CB_SETCURSEL,2,0);
}

static int RefreshFractalType(HWND hWnd, bool refresh = true)
{
	int i = SendDlgItemMessage(hWnd,IDC_COMBO5,CB_GETCURSEL,0,0);
	if (refresh) g_SFT.SetFractalType(i);
	return i;
}

static void UpdateMaxReferences(HWND hWnd)
{
	int i = g_SFT.GetMaxReferences();
	SetDlgItemInt(hWnd,IDC_EDIT3,i,FALSE);
}

static void RefreshMaxReferences(HWND hWnd)
{
	int i = GetDlgItemInt(hWnd,IDC_EDIT3,NULL,FALSE);
	g_SFT.SetMaxReferences(i);
}

static void UpdateGlitchLowTolerance(HWND hWnd)
{
	bool b = g_SFT.GetGlitchLowTolerance();
	SendDlgItemMessage(hWnd,IDC_GLITCHLOWTOLERANCE,BM_SETCHECK,b,0);
}

static void RefreshGlitchLowTolerance(HWND hWnd)
{
	bool b = SendDlgItemMessage(hWnd,IDC_GLITCHLOWTOLERANCE,BM_GETCHECK,0,0);
	g_SFT.SetGlitchLowTolerance(b);
}

static void UpdateApproxLowTolerance(HWND hWnd)
{
	bool b = g_SFT.GetApproxLowTolerance();
	SendDlgItemMessage(hWnd,IDC_CHECK1,BM_SETCHECK,b,0);
}

static void RefreshApproxLowTolerance(HWND hWnd)
{
	bool b = SendDlgItemMessage(hWnd,IDC_CHECK1,BM_GETCHECK,0,0);
	g_SFT.SetApproxLowTolerance(b);
}

static void UpdateAutoApproxTerms(HWND hWnd)
{
	bool b = g_SFT.GetAutoApproxTerms();
	SendDlgItemMessage(hWnd,IDC_CHECK2,BM_SETCHECK,b,0);
}

static void RefreshAutoApproxTerms(HWND hWnd)
{
	bool b = SendDlgItemMessage(hWnd,IDC_CHECK2,BM_GETCHECK,0,0);
	g_SFT.SetAutoApproxTerms(b);
}

static void UpdateApproxTerms(HWND hWnd)
{
	int i = g_SFT.GetApproxTerms();
	SetDlgItemInt(hWnd,IDC_COMBO6,i,FALSE);
}

static void RefreshApproxTerms(HWND hWnd)
{
	char szTerms[256];
	GetDlgItemText(hWnd,IDC_COMBO6,szTerms,sizeof(szTerms));
	int i = atoi(szTerms);
	g_SFT.SetApproxTerms(i);
}

static void UpdateCalculations(HWND hWnd, int nCalc, int nType)
{
  char szCalc[128];
  wsprintf(szCalc,"%d",nCalc);
  int k=strlen(szCalc);
  while(k>3){
    int n=strlen(szCalc);
    while(n>k-4){
      szCalc[n+1]=szCalc[n];
      n--;
    }
    szCalc[k-3]=' ';
    k-=3;
  }
  if(nType)
    strcat(szCalc," 000 000");
  if(g_nPrevCalc!=-1){
    int nC = nCalc-g_nPrevCalc;
    if(!nType)
      nC/=1000000;
    wsprintf(szCalc+strlen(szCalc),", %d M/s",nC);
  }
  g_nPrevCalc=nCalc;
  SetDlgItemText(hWnd,IDC_EDIT8,szCalc);
}

static void UpdateReal(HWND hWnd)
{
	bool b = g_real != 0;
	SendDlgItemMessage(hWnd,IDC_CHECK3,BM_SETCHECK,b,0);
}

static void RefreshReal(HWND hWnd)
{
	g_real = SendDlgItemMessage(hWnd,IDC_CHECK3,BM_GETCHECK,0,0);
}

static void UpdateImag(HWND hWnd)
{
	bool b = g_imag != 0;
	SendDlgItemMessage(hWnd,IDC_CHECK5,BM_SETCHECK,b,0);
}

static void RefreshImag(HWND hWnd)
{
	g_imag = SendDlgItemMessage(hWnd,IDC_CHECK5,BM_GETCHECK,0,0);
}

static void UpdateSeedR(HWND hWnd)
{
	char szTmp[40];
	snprintf(szTmp,40,"%g",g_SeedR);
	SetDlgItemText(hWnd,IDC_EDIT4,szTmp);
}

static void RefreshSeedR(HWND hWnd)
{
	char szTmp[40];
	GetDlgItemText(hWnd,IDC_EDIT4,szTmp,sizeof(szTmp));
	g_SeedR = atof(szTmp);
}

static void UpdateSeedI(HWND hWnd)
{
	char szTmp[40];
	snprintf(szTmp,40,"%g",g_SeedI);
	SetDlgItemText(hWnd,IDC_EDIT9,szTmp);
}

static void RefreshSeedI(HWND hWnd)
{
	char szTmp[40];
	GetDlgItemText(hWnd,IDC_EDIT9,szTmp,sizeof(szTmp));
	g_SeedI = atof(szTmp);
}

static void UpdateFactorAR(HWND hWnd)
{
	char szTmp[40];
	snprintf(szTmp,40,"%g",g_FactorAR);
	SetDlgItemText(hWnd,IDC_EDIT28,szTmp);
}

static void RefreshFactorAR(HWND hWnd)
{
	char szTmp[40];
	GetDlgItemText(hWnd,IDC_EDIT28,szTmp,sizeof(szTmp));
	g_FactorAR = atof(szTmp);
}

static void UpdateFactorAI(HWND hWnd)
{
	char szTmp[40];
	snprintf(szTmp,40,"%g",g_FactorAI);
	SetDlgItemText(hWnd,IDC_EDIT10,szTmp);
}

static void RefreshFactorAI(HWND hWnd)
{
	char szTmp[40];
	GetDlgItemText(hWnd,IDC_EDIT10,szTmp,sizeof(szTmp));
	g_FactorAI = atof(szTmp);
}

extern int WINAPI IterationProc(HWND hWnd,UINT uMsg,WPARAM wParam,LPARAM lParam)
{
	(void) lParam;
	if(uMsg==WM_INITDIALOG || uMsg==WM_TIMER){
		if(uMsg==WM_INITDIALOG){
			SendMessage(hWnd, WM_SETICON, ICON_SMALL, LPARAM(g_hIcon));
			SendMessage(hWnd, WM_SETICON, ICON_BIG, LPARAM(g_hIcon));
			InitToolTip(hWnd,GetModuleHandle(NULL),GetToolText,1);
			UpdateIterations(hWnd);
			SetTimer(hWnd,0,1000,NULL);

			SendDlgItemMessage(hWnd,IDC_COMBO2,CB_ADDSTRING,0,(LPARAM)"High bailout");
			SendDlgItemMessage(hWnd,IDC_COMBO2,CB_ADDSTRING,0,(LPARAM)"Bailout=2");
			UpdateSmoothMethod(hWnd);

			combo5_addstrings(hWnd, IDC_COMBO5);
			SIZE sc;
			HDC hDC = GetDC(NULL);
			SelectObject(hDC,(HFONT)GetStockObject(ANSI_VAR_FONT));
			int i, nMaxWidth=0;
			for(i=0;i<SendDlgItemMessage(hWnd,IDC_COMBO5,CB_GETCOUNT,0,0);i++){
				int n = SendDlgItemMessage(hWnd,IDC_COMBO5,CB_GETLBTEXTLEN,i,0);
				char *szT = new char[n+1];
				SendDlgItemMessage(hWnd,IDC_COMBO5,CB_GETLBTEXT,i,(LPARAM)szT);
				GetTextExtentPoint32A(hDC,szT,strlen(szT),&sc);
				if(sc.cx>nMaxWidth)
					nMaxWidth = sc.cx;
				delete[] szT;
			}
			SendDlgItemMessage(hWnd,IDC_COMBO5,CB_SETDROPPEDWIDTH,nMaxWidth+8+GetSystemMetrics(SM_CXHTHUMB),0);
			UpdateFractalType(hWnd);
      int nType = RefreshFractalType(hWnd, false);

			UpdateGlitchLowTolerance(hWnd);
			UpdateApproxLowTolerance(hWnd);
			UpdateAutoApproxTerms(hWnd);
			UpdatePower(hWnd);
			UpdateMaxReferences(hWnd);

			SendDlgItemMessage(hWnd,IDC_COMBO6,CB_ADDSTRING,0,(LPARAM)"5");
			SendDlgItemMessage(hWnd,IDC_COMBO6,CB_ADDSTRING,0,(LPARAM)"10");
			SendDlgItemMessage(hWnd,IDC_COMBO6,CB_ADDSTRING,0,(LPARAM)"15");
			SendDlgItemMessage(hWnd,IDC_COMBO6,CB_ADDSTRING,0,(LPARAM)"20");
			SendDlgItemMessage(hWnd,IDC_COMBO6,CB_ADDSTRING,0,(LPARAM)"30");
			UpdateApproxTerms(hWnd);

			EnableWindow(GetDlgItem(hWnd,IDC_CHECK2),nType==0 && g_SFT.GetPower()==2);
			EnableWindow(GetDlgItem(hWnd,IDC_COMBO6),nType==0 && g_SFT.GetPower()==2);

			UpdateReal(hWnd);
			UpdateImag(hWnd);
			UpdateSeedR(hWnd);
			UpdateSeedI(hWnd);
			UpdateFactorAR(hWnd);
			UpdateFactorAI(hWnd);
		}
		int nMin, nMax, nCalc=0,nType=0;

		g_SFT.GetIterations(nMin,nMax,&nCalc,&nType);
    UpdateMinIterations(hWnd, nMin);
    UpdateMaxIterations(hWnd, nMax);
		UpdateApproxIterations(hWnd);
		if(uMsg==WM_TIMER){
      UpdateCalculations(hWnd,nCalc,nType);
			if(SendDlgItemMessage(hWnd,IDC_CHECK2,BM_GETCHECK,0,0))
				UpdateApproxTerms(hWnd);
		}
		else
			g_nPrevCalc=-1;
		return 1;
	}
	else if(uMsg==WM_COMMAND){
		if(wParam==IDOK){
			g_bExamineDirty=TRUE;
			RefreshReal(hWnd);
			RefreshImag(hWnd);
			RefreshSeedR(hWnd);
			RefreshSeedI(hWnd);
			RefreshFactorAR(hWnd);
			RefreshFactorAI(hWnd);
			RefreshSmoothMethod(hWnd);
			RefreshApproxTerms(hWnd);
			RefreshGlitchLowTolerance(hWnd);
			RefreshApproxLowTolerance(hWnd);
			RefreshAutoApproxTerms(hWnd);
			RefreshFractalType(hWnd);
			RefreshPower(hWnd);
			RefreshMaxReferences(hWnd);
      RefreshIterations(hWnd);
			ExitToolTip(hWnd);
			EndDialog(hWnd, g_SFT.GetIterations());
		}
		else if(wParam==IDC_CHECK3){
			if(!SendDlgItemMessage(hWnd,IDC_CHECK3,BM_GETCHECK,0,0) && !SendDlgItemMessage(hWnd,IDC_CHECK5,BM_GETCHECK,0,0))
				SendDlgItemMessage(hWnd,IDC_CHECK5,BM_SETCHECK,1,0);
		}
		else if(wParam==IDC_CHECK5){
			if(!SendDlgItemMessage(hWnd,IDC_CHECK3,BM_GETCHECK,0,0) && !SendDlgItemMessage(hWnd,IDC_CHECK5,BM_GETCHECK,0,0))
				SendDlgItemMessage(hWnd,IDC_CHECK3,BM_SETCHECK,1,0);
		}
		else if(wParam==IDCANCEL){
			ExitToolTip(hWnd);
			EndDialog(hWnd,0);
		}
		else if(HIWORD(wParam)==CBN_SELCHANGE && LOWORD(wParam)==IDC_COMBO5){
			int nType = RefreshFractalType(hWnd, false);
			int nPow = RefreshPower(hWnd, false);
			UpdateFractalType(hWnd, nType, nPow);
#ifndef _DEBUG
			EnableWindow(GetDlgItem(hWnd,IDC_COMBO3),nType<=4);
#endif
		}
		int nType = SendDlgItemMessage(hWnd,IDC_COMBO5,CB_GETCURSEL,0,0);
		EnableWindow(GetDlgItem(hWnd,IDC_CHECK2),nType==0 && GetDlgItemInt(hWnd,IDC_COMBO3,NULL,FALSE)==2);
		EnableWindow(GetDlgItem(hWnd,IDC_COMBO6),nType==0 && GetDlgItemInt(hWnd,IDC_COMBO3,NULL,FALSE)==2);
		EnableWindow(GetDlgItem(hWnd,IDC_COMBO6),!SendDlgItemMessage(hWnd,IDC_CHECK2,BM_GETCHECK,0,0));
	}
	return 0;
}

extern const char *IterationToolTip(int nID)
{
  switch(nID){
  case IDC_EDIT1:
    return "Maximum number of iterations";
  case IDC_EDIT2:
    return "Minimum number of iteration in current view";
  case IDC_EDIT5:
    return "Maximum number of iteration in current view";
  case IDC_EDIT7:
    return "Iterations skipped by Series Approximation";
  case IDC_COMBO2:
    return "Bailout value for iterations";
  case IDC_COMBO3:
    return "Power of Mandelbrot function";
  case IDC_EDIT3:
    return "Maximum of extra references for glitch correction";
  case IDC_GLITCHLOWTOLERANCE:
    return "Checked for low tolerance of Glitch Detection\nComplex images may need this to be checked to be rendered correctly\nThe render time may be faster if this checkbox is not checked";
  case IDC_CHECK1:
    return "Checked for low tolerance of Series Approximation\nComplex images may need this to be checked to be rendered correctly\nThe render time may be faster if this checkbox is not checked";
  case IDC_CHECK2:
    return "Terms for Series approximation is adjusted\nbased on the number of pixels to be rendered.";
  case IDC_COMBO5:
    return "List of type of Mandelbrot based Fractals\nSome of them have additional Power options";
  case IDC_COMBO6:
    return "Terms for Series approximation.\nMore terms usually yield more skipped iterations and faster rendering,\nhowever is more time consuming to be processed";
  case IDC_EDIT8:
    return "Display number of calculations performed";
  case IDC_CHECK3:
    return "Include real part when checking bailout.\nUncheck for variation";
  case IDC_CHECK5:
    return "Include imaginary part when checking bailout.\nUncheck for variation";
  case 1002:
    return "Real seed value (0 is standard)";
  case 1003:
    return "Imaginary seed value (0 is standard)";
  case IDOK:
    return "Apply and close";
  case IDCANCEL:
    return "Close and undo";
  default:
    static char tooltip[100];
    snprintf(tooltip, 100, "%d", nID);
    return tooltip;
  }
}