import { NavLink } from "react-router";
import { 
  Activity, 
  Map as MapIcon, 
  Settings, 
  Database,
  FileText,
  BarChart3
} from "lucide-react";

export function Sidebar() {
  const navItems = [
    { name: "Dashboard", path: "/", icon: Activity },
    { name: "Map View", path: "/map", icon: MapIcon },
    { name: "Data Sources", path: "/data", icon: Database },
    { name: "Results", path: "/results", icon: BarChart3 },
    { name: "Reports", path: "/reports", icon: FileText },
    { name: "Settings", path: "/settings", icon: Settings },
  ];

  return (
    <aside className="w-64 bg-[#0f2846] text-slate-300 flex flex-col transition-all hidden md:flex shrink-0">
      <div className="h-14 flex items-center px-4 font-bold text-white tracking-wider border-b border-slate-700/50">
        PFDHA TOOL
      </div>
      <nav className="flex-1 overflow-y-auto py-4 flex flex-col gap-1 px-2">
        {navItems.map((item) => (
          <NavLink
            key={item.name}
            to={item.path}
            className={({ isActive }) =>
              `flex items-center gap-3 px-3 py-2 rounded-md text-sm font-medium transition-colors ${
                isActive
                  ? "bg-blue-600/20 text-white"
                  : "hover:bg-slate-800/50 hover:text-white"
              }`
            }
          >
            <item.icon className="w-4 h-4" />
            {item.name}
          </NavLink>
        ))}
      </nav>
      <div className="p-4 border-t border-slate-700/50 text-xs text-slate-500">
        v2.4.1 (Stable)
      </div>
    </aside>
  );
}
