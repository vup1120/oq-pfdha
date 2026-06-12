import { useState } from "react";
import { ChevronDown, ChevronRight, Play, Settings2 } from "lucide-react";

export function ParameterPanel() {
  const [expandedSection, setExpandedSection] = useState<string>("geometry");

  const toggleSection = (section: string) => {
    setExpandedSection(expandedSection === section ? "" : section);
  };

  return (
    <div className="bg-white rounded-lg border border-slate-200 shadow-sm flex flex-col h-full overflow-hidden">
      <div className="p-4 border-b border-slate-200 bg-slate-50/50 flex items-center gap-2">
        <Settings2 className="w-5 h-5 text-slate-600" />
        <h2 className="font-semibold text-slate-800">Analysis Parameters</h2>
      </div>

      <div className="flex-1 overflow-y-auto p-2">
        {/* Fault Geometry Section */}
        <div className="mb-2">
          <button
            onClick={() => toggleSection("geometry")}
            className="w-full flex items-center justify-between p-2 hover:bg-slate-50 rounded-md transition-colors"
          >
            <span className="font-medium text-sm text-slate-700">Fault Geometry</span>
            {expandedSection === "geometry" ? (
              <ChevronDown className="w-4 h-4 text-slate-400" />
            ) : (
              <ChevronRight className="w-4 h-4 text-slate-400" />
            )}
          </button>
          {expandedSection === "geometry" && (
            <div className="p-2 space-y-4 animate-in fade-in slide-in-from-top-1">
              <div>
                <label className="block text-xs font-medium text-slate-500 mb-1">Fault Type</label>
                <select className="w-full text-sm border border-slate-300 rounded-md p-1.5 bg-white focus:ring-2 focus:ring-blue-500 outline-none">
                  <option>Strike-Slip</option>
                  <option>Normal</option>
                  <option>Reverse</option>
                </select>
              </div>
              <div className="grid grid-cols-2 gap-3">
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Dip (deg)</label>
                  <input type="number" defaultValue="90" className="w-full text-sm border border-slate-300 rounded-md p-1.5 focus:ring-2 focus:ring-blue-500 outline-none" />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Rake (deg)</label>
                  <input type="number" defaultValue="0" className="w-full text-sm border border-slate-300 rounded-md p-1.5 focus:ring-2 focus:ring-blue-500 outline-none" />
                </div>
              </div>
              <div>
                <label className="block text-xs font-medium text-slate-500 mb-1">Seismogenic Depth (km)</label>
                <input type="number" defaultValue="15" className="w-full text-sm border border-slate-300 rounded-md p-1.5 focus:ring-2 focus:ring-blue-500 outline-none" />
              </div>
            </div>
          )}
        </div>

        {/* Site Conditions Section */}
        <div className="mb-2">
          <button
            onClick={() => toggleSection("site")}
            className="w-full flex items-center justify-between p-2 hover:bg-slate-50 rounded-md transition-colors"
          >
            <span className="font-medium text-sm text-slate-700">Site Conditions</span>
            {expandedSection === "site" ? (
              <ChevronDown className="w-4 h-4 text-slate-400" />
            ) : (
              <ChevronRight className="w-4 h-4 text-slate-400" />
            )}
          </button>
          {expandedSection === "site" && (
            <div className="p-2 space-y-4 animate-in fade-in slide-in-from-top-1">
              <div className="grid grid-cols-2 gap-3">
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Lat</label>
                  <input type="number" defaultValue="35.4676" className="w-full text-sm border border-slate-300 rounded-md p-1.5 focus:ring-2 focus:ring-blue-500 outline-none" />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Lon</label>
                  <input type="number" defaultValue="-97.5164" className="w-full text-sm border border-slate-300 rounded-md p-1.5 focus:ring-2 focus:ring-blue-500 outline-none" />
                </div>
              </div>
              <div>
                <label className="block text-xs font-medium text-slate-500 mb-1">Distance to Fault (km)</label>
                <input type="number" defaultValue="2.5" className="w-full text-sm border border-slate-300 rounded-md p-1.5 focus:ring-2 focus:ring-blue-500 outline-none" />
              </div>
            </div>
          )}
        </div>

        {/* Models Section */}
        <div className="mb-2">
          <button
            onClick={() => toggleSection("models")}
            className="w-full flex items-center justify-between p-2 hover:bg-slate-50 rounded-md transition-colors"
          >
            <span className="font-medium text-sm text-slate-700">Displacement Models</span>
            {expandedSection === "models" ? (
              <ChevronDown className="w-4 h-4 text-slate-400" />
            ) : (
              <ChevronRight className="w-4 h-4 text-slate-400" />
            )}
          </button>
          {expandedSection === "models" && (
            <div className="p-2 space-y-3 animate-in fade-in slide-in-from-top-1">
              <label className="flex items-center gap-2 cursor-pointer">
                <input type="checkbox" defaultChecked className="rounded text-blue-600 focus:ring-blue-500 w-4 h-4" />
                <span className="text-sm text-slate-700">Petersen et al. (2011)</span>
              </label>
              <label className="flex items-center gap-2 cursor-pointer">
                <input type="checkbox" defaultChecked className="rounded text-blue-600 focus:ring-blue-500 w-4 h-4" />
                <span className="text-sm text-slate-700">Youngs et al. (2003)</span>
              </label>
              <label className="flex items-center gap-2 cursor-pointer">
                <input type="checkbox" className="rounded text-blue-600 focus:ring-blue-500 w-4 h-4" />
                <span className="text-sm text-slate-700">Moss et al. (2015)</span>
              </label>
            </div>
          )}
        </div>
      </div>

      <div className="p-4 border-t border-slate-200 bg-slate-50/50">
        <button className="w-full bg-[#0f2846] hover:bg-[#1a3d66] text-white flex items-center justify-center gap-2 py-2.5 rounded-md font-medium text-sm transition-colors shadow-sm">
          <Play className="w-4 h-4" />
          Run Analysis
        </button>
      </div>
    </div>
  );
}
