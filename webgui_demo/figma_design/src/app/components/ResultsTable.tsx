import { Download, Filter, MoreHorizontal } from "lucide-react";

const tableData = [
  { returnPeriod: "475 yrs", exceedance: "2.11e-3", petersenDisp: "0.85", youngsDisp: "0.92", meanDisp: "0.88" },
  { returnPeriod: "975 yrs", exceedance: "1.03e-3", petersenDisp: "1.24", youngsDisp: "1.45", meanDisp: "1.34" },
  { returnPeriod: "2,475 yrs", exceedance: "4.04e-4", petersenDisp: "2.56", youngsDisp: "3.10", meanDisp: "2.83" },
  { returnPeriod: "4,975 yrs", exceedance: "2.01e-4", petersenDisp: "4.12", youngsDisp: "4.85", meanDisp: "4.48" },
  { returnPeriod: "9,975 yrs", exceedance: "1.00e-4", petersenDisp: "6.35", youngsDisp: "7.20", meanDisp: "6.77" },
];

export function ResultsTable() {
  return (
    <div className="bg-white rounded-lg border border-slate-200 shadow-sm flex flex-col h-full overflow-hidden">
      <div className="p-4 border-b border-slate-200 flex items-center justify-between bg-slate-50/50">
        <h3 className="font-semibold text-slate-800 text-sm">Displacement Summary</h3>
        <div className="flex items-center gap-2">
          <button className="flex items-center gap-1.5 px-3 py-1.5 border border-slate-300 rounded-md text-xs font-medium text-slate-700 hover:bg-slate-50 transition-colors">
            <Filter className="w-3.5 h-3.5" />
            Filter
          </button>
          <button className="flex items-center gap-1.5 px-3 py-1.5 border border-slate-300 rounded-md text-xs font-medium text-slate-700 hover:bg-slate-50 transition-colors">
            <Download className="w-3.5 h-3.5" />
            Export CSV
          </button>
        </div>
      </div>
      
      <div className="flex-1 overflow-auto">
        <table className="w-full text-sm text-left">
          <thead className="text-xs text-slate-600 bg-slate-50 sticky top-0 z-10">
            <tr>
              <th className="px-4 py-3 font-medium border-b border-slate-200">Return Period</th>
              <th className="px-4 py-3 font-medium border-b border-slate-200">Exceedance Rate</th>
              <th className="px-4 py-3 font-medium border-b border-slate-200 text-right">Petersen Disp. (cm)</th>
              <th className="px-4 py-3 font-medium border-b border-slate-200 text-right">Youngs Disp. (cm)</th>
              <th className="px-4 py-3 font-medium border-b border-slate-200 text-right">Mean Disp. (cm)</th>
              <th className="px-4 py-3 font-medium border-b border-slate-200 text-center w-10"></th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-100">
            {tableData.map((row, idx) => (
              <tr key={idx} className="hover:bg-slate-50/80 transition-colors">
                <td className="px-4 py-3 font-medium text-slate-800">{row.returnPeriod}</td>
                <td className="px-4 py-3 text-slate-600 font-mono text-xs">{row.exceedance}</td>
                <td className="px-4 py-3 text-right text-slate-700">{row.petersenDisp}</td>
                <td className="px-4 py-3 text-right text-slate-700">{row.youngsDisp}</td>
                <td className="px-4 py-3 text-right font-medium text-[#0f2846]">{row.meanDisp}</td>
                <td className="px-4 py-3 text-center">
                  <button className="text-slate-400 hover:text-slate-600">
                    <MoreHorizontal className="w-4 h-4" />
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
