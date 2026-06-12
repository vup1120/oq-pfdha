import { createBrowserRouter } from "react-router";
import { RootLayout } from "./components/RootLayout";
import { Dashboard } from "./components/Dashboard";

export const router = createBrowserRouter([
  {
    path: "/",
    Component: RootLayout,
    children: [
      {
        index: true,
        Component: Dashboard,
      },
      // Mock other routes
      {
        path: "*",
        Component: () => <div className="p-8 text-slate-500">View not implemented in this demo.</div>
      }
    ],
  },
]);
