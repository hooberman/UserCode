{
//=========Macro generated from canvas: c/the canvas
//=========  (Fri Feb  4 12:59:11 2011) by ROOT version5.27/06
   TCanvas *c = new TCanvas("c", "the canvas",5,124,650,500);
   gStyle->SetOptStat(0);
   c->Range(-0.128266,0.07692306,1.059382,1.102564);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#f0f0f0");
   c->SetFillColor(ci);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetGridx();
   c->SetGridy();
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetRightMargin(0.05);
   c->SetBottomMargin(0.12);

   ci = TColor::GetColor("#fffffd");
   c->SetFrameFillColor(ci);
   c->SetFrameBorderMode(0);

   ci = TColor::GetColor("#fffffd");
   c->SetFrameFillColor(ci);
   c->SetFrameBorderMode(0);
   
   TH2F *frame = new TH2F("frame","Background rejection versus Signal efficiency",500,0,1,500,0.2,1);
   frame->SetStats(0);
   frame->SetLineWidth(2);
   frame->SetMarkerStyle(21);
   frame->SetMarkerSize(0.3);
   frame->GetXaxis()->SetTitle("Signal efficiency");
   frame->GetXaxis()->SetLabelOffset(0.012);
   frame->GetXaxis()->SetTitleSize(0.045);
   frame->GetXaxis()->SetTitleOffset(1.25);
   frame->GetYaxis()->SetTitle("Background rejection");
   frame->GetYaxis()->SetLabelOffset(0.012);
   frame->GetYaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleOffset(1.22);
   frame->Draw("");
   
   TH1F *MVA_Fisher_rejBvsS = new TH1F("MVA_Fisher_rejBvsS","MVA_Fisher",100,0,1);
   MVA_Fisher_rejBvsS->SetBinContent(1,0.9985881);
   MVA_Fisher_rejBvsS->SetBinContent(2,0.9961708);
   MVA_Fisher_rejBvsS->SetBinContent(3,0.9940513);
   MVA_Fisher_rejBvsS->SetBinContent(4,0.9910768);
   MVA_Fisher_rejBvsS->SetBinContent(5,0.9881001);
   MVA_Fisher_rejBvsS->SetBinContent(6,0.9870523);
   MVA_Fisher_rejBvsS->SetBinContent(7,0.9830642);
   MVA_Fisher_rejBvsS->SetBinContent(8,0.9822509);
   MVA_Fisher_rejBvsS->SetBinContent(9,0.9803713);
   MVA_Fisher_rejBvsS->SetBinContent(10,0.9786204);
   MVA_Fisher_rejBvsS->SetBinContent(11,0.9766313);
   MVA_Fisher_rejBvsS->SetBinContent(12,0.9735237);
   MVA_Fisher_rejBvsS->SetBinContent(13,0.9704745);
   MVA_Fisher_rejBvsS->SetBinContent(14,0.968612);
   MVA_Fisher_rejBvsS->SetBinContent(15,0.9660085);
   MVA_Fisher_rejBvsS->SetBinContent(16,0.9638223);
   MVA_Fisher_rejBvsS->SetBinContent(17,0.960419);
   MVA_Fisher_rejBvsS->SetBinContent(18,0.957853);
   MVA_Fisher_rejBvsS->SetBinContent(19,0.955829);
   MVA_Fisher_rejBvsS->SetBinContent(20,0.9542154);
   MVA_Fisher_rejBvsS->SetBinContent(21,0.9493043);
   MVA_Fisher_rejBvsS->SetBinContent(22,0.9463041);
   MVA_Fisher_rejBvsS->SetBinContent(23,0.9431273);
   MVA_Fisher_rejBvsS->SetBinContent(24,0.9400573);
   MVA_Fisher_rejBvsS->SetBinContent(25,0.9370363);
   MVA_Fisher_rejBvsS->SetBinContent(26,0.9344524);
   MVA_Fisher_rejBvsS->SetBinContent(27,0.9304616);
   MVA_Fisher_rejBvsS->SetBinContent(28,0.9276335);
   MVA_Fisher_rejBvsS->SetBinContent(29,0.9243718);
   MVA_Fisher_rejBvsS->SetBinContent(30,0.9218449);
   MVA_Fisher_rejBvsS->SetBinContent(31,0.9187946);
   MVA_Fisher_rejBvsS->SetBinContent(32,0.9155319);
   MVA_Fisher_rejBvsS->SetBinContent(33,0.9124193);
   MVA_Fisher_rejBvsS->SetBinContent(34,0.9084308);
   MVA_Fisher_rejBvsS->SetBinContent(35,0.9060782);
   MVA_Fisher_rejBvsS->SetBinContent(36,0.90308);
   MVA_Fisher_rejBvsS->SetBinContent(37,0.8986884);
   MVA_Fisher_rejBvsS->SetBinContent(38,0.8943672);
   MVA_Fisher_rejBvsS->SetBinContent(39,0.8914782);
   MVA_Fisher_rejBvsS->SetBinContent(40,0.888705);
   MVA_Fisher_rejBvsS->SetBinContent(41,0.8856748);
   MVA_Fisher_rejBvsS->SetBinContent(42,0.8803468);
   MVA_Fisher_rejBvsS->SetBinContent(43,0.8774128);
   MVA_Fisher_rejBvsS->SetBinContent(44,0.8750042);
   MVA_Fisher_rejBvsS->SetBinContent(45,0.8686271);
   MVA_Fisher_rejBvsS->SetBinContent(46,0.8631946);
   MVA_Fisher_rejBvsS->SetBinContent(47,0.8592953);
   MVA_Fisher_rejBvsS->SetBinContent(48,0.854953);
   MVA_Fisher_rejBvsS->SetBinContent(49,0.8523598);
   MVA_Fisher_rejBvsS->SetBinContent(50,0.8477351);
   MVA_Fisher_rejBvsS->SetBinContent(51,0.8433961);
   MVA_Fisher_rejBvsS->SetBinContent(52,0.8383685);
   MVA_Fisher_rejBvsS->SetBinContent(53,0.8329851);
   MVA_Fisher_rejBvsS->SetBinContent(54,0.828046);
   MVA_Fisher_rejBvsS->SetBinContent(55,0.8220918);
   MVA_Fisher_rejBvsS->SetBinContent(56,0.8179269);
   MVA_Fisher_rejBvsS->SetBinContent(57,0.814643);
   MVA_Fisher_rejBvsS->SetBinContent(58,0.8106885);
   MVA_Fisher_rejBvsS->SetBinContent(59,0.8077152);
   MVA_Fisher_rejBvsS->SetBinContent(60,0.8024303);
   MVA_Fisher_rejBvsS->SetBinContent(61,0.796521);
   MVA_Fisher_rejBvsS->SetBinContent(62,0.7902446);
   MVA_Fisher_rejBvsS->SetBinContent(63,0.7864058);
   MVA_Fisher_rejBvsS->SetBinContent(64,0.7815406);
   MVA_Fisher_rejBvsS->SetBinContent(65,0.7756431);
   MVA_Fisher_rejBvsS->SetBinContent(66,0.7691262);
   MVA_Fisher_rejBvsS->SetBinContent(67,0.7636815);
   MVA_Fisher_rejBvsS->SetBinContent(68,0.7578233);
   MVA_Fisher_rejBvsS->SetBinContent(69,0.7511776);
   MVA_Fisher_rejBvsS->SetBinContent(70,0.7436728);
   MVA_Fisher_rejBvsS->SetBinContent(71,0.7363932);
   MVA_Fisher_rejBvsS->SetBinContent(72,0.7290592);
   MVA_Fisher_rejBvsS->SetBinContent(73,0.72113);
   MVA_Fisher_rejBvsS->SetBinContent(74,0.7134986);
   MVA_Fisher_rejBvsS->SetBinContent(75,0.7064184);
   MVA_Fisher_rejBvsS->SetBinContent(76,0.7007454);
   MVA_Fisher_rejBvsS->SetBinContent(77,0.693757);
   MVA_Fisher_rejBvsS->SetBinContent(78,0.6856257);
   MVA_Fisher_rejBvsS->SetBinContent(79,0.6796134);
   MVA_Fisher_rejBvsS->SetBinContent(80,0.6723097);
   MVA_Fisher_rejBvsS->SetBinContent(81,0.6652589);
   MVA_Fisher_rejBvsS->SetBinContent(82,0.6565596);
   MVA_Fisher_rejBvsS->SetBinContent(83,0.6481721);
   MVA_Fisher_rejBvsS->SetBinContent(84,0.6374931);
   MVA_Fisher_rejBvsS->SetBinContent(85,0.6284167);
   MVA_Fisher_rejBvsS->SetBinContent(86,0.6192074);
   MVA_Fisher_rejBvsS->SetBinContent(87,0.6061083);
   MVA_Fisher_rejBvsS->SetBinContent(88,0.598141);
   MVA_Fisher_rejBvsS->SetBinContent(89,0.5899556);
   MVA_Fisher_rejBvsS->SetBinContent(90,0.5774851);
   MVA_Fisher_rejBvsS->SetBinContent(91,0.566143);
   MVA_Fisher_rejBvsS->SetBinContent(92,0.5560098);
   MVA_Fisher_rejBvsS->SetBinContent(93,0.5431424);
   MVA_Fisher_rejBvsS->SetBinContent(94,0.530799);
   MVA_Fisher_rejBvsS->SetBinContent(95,0.513099);
   MVA_Fisher_rejBvsS->SetBinContent(96,0.4930019);
   MVA_Fisher_rejBvsS->SetBinContent(97,0.4737727);
   MVA_Fisher_rejBvsS->SetBinContent(98,0.4512871);
   MVA_Fisher_rejBvsS->SetBinContent(99,0.4288951);
   MVA_Fisher_rejBvsS->SetBinContent(100,0.387387);
   MVA_Fisher_rejBvsS->SetEntries(100);
   MVA_Fisher_rejBvsS->SetLineWidth(3);
   MVA_Fisher_rejBvsS->GetXaxis()->SetTitle("Signal eff");
   MVA_Fisher_rejBvsS->GetXaxis()->SetLabelSize(0.05);
   MVA_Fisher_rejBvsS->GetXaxis()->SetTitleSize(0.055);
   MVA_Fisher_rejBvsS->GetYaxis()->SetTitle("Backgr rejection (1-eff)");
   MVA_Fisher_rejBvsS->GetYaxis()->SetLabelSize(0.05);
   MVA_Fisher_rejBvsS->GetYaxis()->SetTitleSize(0.055);
   MVA_Fisher_rejBvsS->GetZaxis()->SetLabelSize(0.05);
   MVA_Fisher_rejBvsS->GetZaxis()->SetTitleSize(0.055);
   MVA_Fisher_rejBvsS->Draw("csame");
   
   TH2F *frame = new TH2F("frame","Background rejection versus Signal efficiency",500,0,1,500,0.2,1);
   frame->SetStats(0);
   frame->SetLineWidth(2);
   frame->SetMarkerStyle(21);
   frame->SetMarkerSize(0.3);
   frame->GetXaxis()->SetTitle("Signal efficiency");
   frame->GetXaxis()->SetLabelOffset(0.012);
   frame->GetXaxis()->SetTitleSize(0.045);
   frame->GetXaxis()->SetTitleOffset(1.25);
   frame->GetYaxis()->SetTitle("Background rejection");
   frame->GetYaxis()->SetLabelOffset(0.012);
   frame->GetYaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleOffset(1.22);
   frame->Draw("sameaxis");
   
   TLegend *leg = new TLegend(0.15,0.171,0.5,0.281,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);

   ci = TColor::GetColor("#7d8b9d");
   leg->SetLineColor(ci);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","MVA Method:","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("MVA_Fisher_rejBvsS","Fisher","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.9390678,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#5d6b7d");
   pt->SetFillColor(ci);

   ci = TColor::GetColor("#7d8b9d");
   pt->SetLineColor(ci);
   pt->SetTextColor(0);
   TText *text = pt->AddText("Background rejection versus Signal efficiency");
   pt->Draw();
  
// ------------>Primitives in pad: imgpad
   TPad *imgpad = new TPad("imgpad", "imgpad",0.789257,0.91,0.95,0.965);
   imgpad->Draw();
   imgpad->cd();
   imgpad->Range(0,0,1,1);
   imgpad->SetFillColor(0);
   imgpad->SetBorderMode(0);
   imgpad->SetBorderSize(2);
   imgpad->SetTickx(1);
   imgpad->SetTicky(1);
   imgpad->SetLeftMargin(0);
   imgpad->SetRightMargin(0);
   imgpad->SetTopMargin(0);
   imgpad->SetBottomMargin(0);

   ci = TColor::GetColor("#fffffd");
   imgpad->SetFrameFillColor(ci);
   imgpad->SetFrameBorderMode(0);
   imgpad->SetFrameLineColor(0);
   imgpad->SetFrameBorderMode(0);

/* XPM */
 char *xpm_tmva_logo_gif_1[] = {
/* columns rows colors chars-per-pixel */
"103 25 89 1",
"  c #6A4F7B",
". c #67557B",
"X c #75597A",
"o c #FF0303",
"O c #FF1919",
"+ c #FF2A09",
"@ c #FF3D0C",
"# c #FF2020",
"$ c #FF3A3A",
"% c #FF6012",
"& c #FF7B18",
"* c #B17B58",
"= c #8B646A",
"- c #936964",
"; c #FF4444",
": c #FF5252",
"> c #FF6D6D",
", c #FF7E6A",
"< c #FF7979",
"1 c #34399B",
"2 c #3536A0",
"3 c #423C91",
"4 c #3F409E",
"5 c #49408D",
"6 c #5E4B81",
"7 c #464197",
"8 c #624E81",
"9 c #695180",
"0 c #4444A5",
"q c #5D5DAD",
"w c #5D5DB3",
"e c #6C62A3",
"r c #6D6DB9",
"t c #7474B8",
"y c #7D7DC3",
"u c #D37C95",
"i c #FF881B",
"p c #FF9E1F",
"a c #D28E3F",
"s c #FF8721",
"d c #F59F2A",
"f c #E79833",
"g c #FFA621",
"h c #FFAE30",
"j c #FFB237",
"k c #BE824C",
"l c #CD8B45",
"z c #D38D40",
"x c #D79540",
"c c #FFB845",
"v c #FFBE54",
"b c #F2BC6D",
"n c #FFC25D",
"m c #FFC360",
"M c #FFC673",
"N c #9183AB",
"B c #8B81B2",
"V c #FF8686",
"C c #FFAFAF",
"Z c #FFBFBF",
"A c #8686C6",
"S c #9C9CCC",
"D c #9B9BD1",
"F c #A0A0D3",
"G c #B9B9DB",
"H c #BEBEE2",
"J c #C6BBC9",
"K c #D3B2CC",
"L c #FFD48B",
"P c #FFDB9E",
"I c #E3CBBF",
"U c #FFDEAF",
"Y c #FFDFB7",
"T c #FFE1AC",
"R c #FFE8BD",
"E c #D3C2C1",
"W c #FFC2C2",
"Q c #F9D9DF",
"! c #CCCCEB",
"~ c #D3CCE5",
"^ c #D0D0EC",
"/ c #DDDDF1",
"( c #FFEBC3",
") c #FFEDD2",
"_ c #FFF0D3",
"` c #E7E7F6",
"' c #FFF8E4",
"] c #FFFFFF",
"[ c None",
"]]]]]]]]]]]]]]]]]](mggjM_]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
"]]]]]]]]]]]]]]]]]ngggggggT]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
"]]]]]]]]]]]]]]]]LgggggggggR]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
"]]]]]]]]]]]]]]]'ggggggggggc]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
"]]]]]]]]]]]]]]]Lggggggggggg_]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
"]]]]]]]]]]]]]]]ngggggggggggP]]]]]]]]]]]]]]WWZWWWZWWWZQ^!^^]]]]]]^!^^]]!^!]]]]]]]!^!]]]]]/^!^]]]]]]]]]]]",
"]]]]]]]]]]]]]]]ngggggggggggP]]]]]]]]]]]]]Wooooooooooou1124]]]]]F112q]`112S]]]]]S12t]]]]]0112F]]]]]]]]]]",
"]]]]]]]]]]]]]]]Lggggggggggg(]]]]]]]]]]]]]<oooooooooooK1211F]]]]2112q]]q210]]]]]q12G]]]]H1121q]]]]]]]]]]",
"]]]]]]]]]]]]]]](ggggggggggg]]]]]]]]]]]]]]QWWZ<oooCWWW~11210]]]G1211w]]G211!]]]!110]]]]]r11211/]]]]]]]]]",
"]]]]]]]]]]]]]]]]vgggggggggP]]]]]]]]]]]]]]]]]]>oo#]]]]^12111F]]r1121w]]]211D]]]A11y]]]]]121011D]]]]]]]]]",
"]]]]]]]]]]]]]]]]'jgggggggv]]]]]]]]]]]]]]]]]]];oo>]]]]^121q1w]!11w11w]]]t12q]]]011H]]]]F11yG11w]]]]]]]]]",
"]]]]]]]]]]]]]]]](MLggggjUL']]]]]]]]]]]]]]]]]]oooC]]]]^121S21!q10A11w]]]H112`]G210]]]]]r11!]121^]]]]]]]]",
"]]]]]]]]'((((_]]M]]'(U(]]]m]](((((]]]]]]]]]]WoooQ]]]]^121`11012GA11w]]]]121F]t21A]]]]]210]]r11A]]]]]]]]",
"]]]]]]]_']]]]](M]]]]]]]]]]'v(]]]]](']]]]]]]]VooO]]]]]^121]t2110]A11w]]]]r11w]121!]]]]F11t]]D110]]]]]]]]",
"]]]]]]_]]]]]]]](]]]]]]]]]]]_]]]]]]]'']]]]]]]$oo:]]]]]^121]/111A]A11w]]]]H112S110]]]]]q210qr0112!]]]]]]]",
")((((R]]]]]]]]]](((((((((((]]]]]]]]]((((((((ooo,(((((E111((e11I(B11q(((((711211N((((J1121112111t((((((_",
"cggggL]]]]]]]]]]Uggggggggg']]]]]]]]]]cgggggiooo&gggggz111ggggggg=126ggggg-11211lgggg=123.9 96125ggggggv",
"cgggg(]]]]]]]]]]_ggggggggc]]]]]]]]]]]Lggggg%ooopgggggz111ggggggg=126gggggf11217ggggd111*ggggf112lgggggv",
"cgggg(]]]]]]]]]]_ggggggggv]]]]]]]]]]]Pggggg@oo+ggggggz111ggggggg=126gggggg8111-ggggk123gggggg611Xgggggv",
"cgggg(]]]]]]]]]]_ggggggggj]]]]]]]]]]]Mgggggi&isgggggggzlagggggggflaxggggggflalgggggflaxggggggflalgggggn",
"cggggP]]]]]]]]]]Rggggggggg']]]]]]]]]]cggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggv",
"cggggj]]]]]]]]]]ngggggggggP]]]]]]]]](gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggn",
"cgggggv]]]]]]]]Pggggggggggg_]]]]]]]'hgggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggn",
"cggggggm]]]]]]Lggggggggggggg(]]]]](cggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggv",
"bnvnvnvnmY((RMnvnvnvnvnvnvnvnP(((UnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvnvL"
};


   TImage *tmva_logo_gif_ = TImage::Create();
   tmva_logo_gif_->SetImageBuffer(xpm_tmva_logo_gif_1, TImage::kXpm);
   tmva_logo_gif_->Draw();
   imgpad->Modified();
   c->cd();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
