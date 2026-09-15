import AppKit
import Foundation
@preconcurrency import WebKit

final class AppDelegate: NSObject, NSApplicationDelegate, NSWindowDelegate, WKNavigationDelegate, WKUIDelegate, WKDownloadDelegate {
    private var child: Process?
    private var window: NSWindow!
    private var webView: WKWebView!
    private var startupTimer: Timer?
    private var loadedAddress: String?
    private var isQuitting = false
    private var leaveServiceRunning = false
    private var terminationPending = false
    private var downloads: [ObjectIdentifier: WKDownload] = [:]
    private var downloadFiles: [ObjectIdentifier: URL] = [:]
    private var savePanels = 0
    private let support: URL = {
        let arguments = CommandLine.arguments
        if let index = arguments.firstIndex(of: "--home"), index + 1 < arguments.count {
            return URL(fileURLWithPath: arguments[index + 1], isDirectory: true).standardizedFileURL
        }
        return FileManager.default.homeDirectoryForCurrentUser.appendingPathComponent("Library/Application Support/DynaMol")
    }()

    func applicationDidFinishLaunching(_ notification: Notification) {
        makeMenus()
        let configuration = WKWebViewConfiguration()
        configuration.preferences.isElementFullscreenEnabled = true
        webView = WKWebView(frame: .zero, configuration: configuration)
        webView.navigationDelegate = self
        webView.uiDelegate = self
        webView.allowsBackForwardNavigationGestures = false
        window = NSWindow(contentRect: NSRect(x: 0, y: 0, width: 1360, height: 900),
                          styleMask: [.titled, .closable, .miniaturizable, .resizable], backing: .buffered, defer: false)
        window.title = "DynaMol"
        window.minSize = NSSize(width: 1000, height: 700)
        window.contentView = webView
        window.delegate = self
        window.isReleasedWhenClosed = false
        window.collectionBehavior = [.fullScreenPrimary]
        window.setFrameAutosaveName("DynaMolWorkspace")
        window.center()
        showStartingPage()
        openWorkspace()
        start()
    }

    private func menuItem(_ title: String, _ action: Selector, _ key: String = "", target: AnyObject? = nil) -> NSMenuItem {
        let result = NSMenuItem(title: title, action: action, keyEquivalent: key)
        result.target = target
        return result
    }

    private func makeMenus() {
        let main = NSMenu()
        func addMenu(_ title: String) -> NSMenu {
            let root = NSMenuItem(title: title, action: nil, keyEquivalent: "")
            let menu = NSMenu(title: title)
            root.submenu = menu
            main.addItem(root)
            return menu
        }
        let appMenu = addMenu("DynaMol")
        appMenu.addItem(menuItem("About DynaMol", #selector(NSApplication.orderFrontStandardAboutPanel(_:))))
        appMenu.addItem(.separator())
        let services = NSMenu(title: "Services")
        let servicesItem = NSMenuItem(title: "Services", action: nil, keyEquivalent: "")
        servicesItem.submenu = services
        appMenu.addItem(servicesItem)
        NSApp.servicesMenu = services
        appMenu.addItem(.separator())
        appMenu.addItem(menuItem("Hide DynaMol", #selector(NSApplication.hide(_:)), "h"))
        let hideOthers = menuItem("Hide Others", #selector(NSApplication.hideOtherApplications(_:)), "h")
        hideOthers.keyEquivalentModifierMask = [.command, .option]
        appMenu.addItem(hideOthers)
        appMenu.addItem(menuItem("Show All", #selector(NSApplication.unhideAllApplications(_:))))
        appMenu.addItem(.separator())
        appMenu.addItem(menuItem("Quit DynaMol", #selector(NSApplication.terminate(_:)), "q"))

        let file = addMenu("File")
        file.addItem(menuItem("Open DynaMol Window", #selector(openWorkspace), "o", target: self))
        file.addItem(menuItem("Show Workspace Files", #selector(showFiles), target: self))
        file.addItem(.separator())
        file.addItem(menuItem("Close Window", #selector(NSWindow.performClose(_:)), "w"))
        let edit = addMenu("Edit")
        edit.addItem(menuItem("Undo", Selector(("undo:")), "z"))
        let redo = menuItem("Redo", Selector(("redo:")), "z")
        redo.keyEquivalentModifierMask = [.command, .shift]
        edit.addItem(redo)
        edit.addItem(.separator())
        edit.addItem(menuItem("Cut", #selector(NSText.cut(_:)), "x"))
        edit.addItem(menuItem("Copy", #selector(NSText.copy(_:)), "c"))
        edit.addItem(menuItem("Paste", #selector(NSText.paste(_:)), "v"))
        edit.addItem(menuItem("Select All", #selector(NSText.selectAll(_:)), "a"))
        let view = addMenu("View")
        view.addItem(menuItem("Reload Workspace", #selector(reloadWorkspace), "r", target: self))
        let fullscreen = menuItem("Enter Full Screen", #selector(NSWindow.toggleFullScreen(_:)), "f")
        fullscreen.keyEquivalentModifierMask = [.command, .control]
        view.addItem(fullscreen)
        let windows = addMenu("Window")
        windows.addItem(menuItem("Minimize", #selector(NSWindow.performMiniaturize(_:)), "m"))
        windows.addItem(menuItem("Zoom", #selector(NSWindow.performZoom(_:))))
        windows.addItem(.separator())
        windows.addItem(menuItem("Bring All to Front", #selector(NSApplication.arrangeInFront(_:))))
        NSApp.windowsMenu = windows
        let help = addMenu("Help")
        help.addItem(menuItem("Show DynaMol Logs", #selector(showLogs), target: self))
        NSApp.helpMenu = help
        NSApp.mainMenu = main
    }

    private func showStartingPage() {
        loadedAddress = nil
        webView.loadHTMLString("""
        <!doctype html><html><meta name="viewport" content="width=device-width,initial-scale=1"><style>
        body{margin:0;background:#0b121b;color:#e3f0ed;font:16px -apple-system,BlinkMacSystemFont,sans-serif;display:grid;place-items:center;min-height:100vh}.card{width:min(460px,80vw);padding:42px;border:1px solid #294139;border-radius:22px;background:linear-gradient(130deg,#142a29,#121e2a)}h1{font-size:34px;letter-spacing:-1px}p{color:#b7c8cc;line-height:1.65}.wheel{border:3px solid #28483f;border-top-color:#83e5c5;width:28px;height:28px;border-radius:50%;animation:spin 1s linear infinite}@keyframes spin{to{transform:rotate(360deg)}}small{color:#91aaa4}</style><main class="card"><div class="wheel"></div><h1>Dyna<span style="color:#83e5c5">Mol</span></h1><p>Opening your molecular workspace…</p><small>OpenMM and GROMACS are included. Your molecules stay on this Mac.</small></main></html>
        """, baseURL: nil)
    }

    private func start() {
        guard child?.isRunning != true, !isQuitting, let resources = Bundle.main.resourceURL else { return }
        let process = Process()
        process.executableURL = resources.appendingPathComponent("python/bin/python3.12")
        var arguments = Array(CommandLine.arguments.dropFirst())
        if !arguments.contains("--no-browser") { arguments.append("--no-browser") }
        process.arguments = ["-I", "-B", resources.appendingPathComponent("bootstrap.py").path] + arguments
        let logs = support.appendingPathComponent("Logs")
        try? FileManager.default.createDirectory(at: logs, withIntermediateDirectories: true)
        let file = logs.appendingPathComponent("desktop.log")
        if !FileManager.default.fileExists(atPath: file.path) { FileManager.default.createFile(atPath: file.path, contents: nil) }
        if let handle = try? FileHandle(forWritingTo: file) {
            _ = try? handle.seekToEnd()
            process.standardOutput = handle
            process.standardError = handle
        }
        process.terminationHandler = { [weak self] process in
            DispatchQueue.main.async {
                guard let self = self else { return }
                if self.child === process { self.child = nil }
                guard !self.isQuitting else { return }
                // A second launcher exits normally when it reuses an existing service.
                if process.terminationStatus != 0 {
                    self.startupTimer?.invalidate()
                    self.showError("DynaMol could not open", "Your saved molecules are safe. Details are in \(logs.path). You can reopen the workspace with File → Open DynaMol Window.")
                } else { self.followStartup() }
            }
        }
        do {
            try process.run()
            child = process
            startupTimer?.invalidate()
            startupTimer = Timer.scheduledTimer(withTimeInterval: 0.4, repeats: true) { [weak self] _ in self?.followStartup() }
            followStartup()
        } catch { showError("DynaMol could not start", error.localizedDescription) }
    }

    private func stateFile(_ name: String) -> [String: Any]? {
        guard let data = try? Data(contentsOf: support.appendingPathComponent(name)),
              let object = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
              let address = object["url"] as? String, let url = URL(string: address),
              url.scheme == "http", url.host == "127.0.0.1", url.port != nil,
              let pid = object["pid"] as? Int, pid > 0, kill(pid_t(pid), 0) == 0 else { return nil }
        return object
    }

    private func serviceState() -> [String: Any]? { stateFile("service.json") }

    private func followStartup() {
        guard !isQuitting else { return }
        let service = serviceState()
        if let object = service ?? stateFile("startup.json"),
           let address = object["url"] as? String, let url = URL(string: address) {
            if loadedAddress != address {
                loadedAddress = address
                webView.load(URLRequest(url: url))
            }
            if service != nil { startupTimer?.invalidate(); startupTimer = nil }
        }
    }

    @objc func openWorkspace() {
        window?.makeKeyAndOrderFront(nil)
        NSApp.activate(ignoringOtherApps: true)
        if webView != nil {
            followStartup()
            if child == nil && serviceState() == nil && stateFile("startup.json") == nil { start() }
        }
    }
    @objc func reloadWorkspace() {
        if serviceState() == nil && child == nil { showStartingPage(); start() }
        else { webView.reload() }
    }
    @objc func showFiles() { NSWorkspace.shared.open(support.appendingPathComponent("Workspace")) }
    @objc func showLogs() { NSWorkspace.shared.open(support.appendingPathComponent("Logs")) }
    func applicationShouldHandleReopen(_ sender: NSApplication, hasVisibleWindows flag: Bool) -> Bool { openWorkspace(); return true }
    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { false }

    private func activeJobs(_ completion: @escaping ([[String: Any]]?) -> Void) {
        guard let state = serviceState(), let address = state["url"] as? String,
              let base = URL(string: address) else { completion(nil); return }
        var request = URLRequest(url: base.appendingPathComponent("api/jobs")); request.timeoutInterval = 3
        URLSession.shared.dataTask(with: request) { data, response, _ in
            guard (response as? HTTPURLResponse)?.statusCode == 200, let data = data,
                  let jobs = try? JSONSerialization.jsonObject(with: data) as? [[String: Any]] else {
                DispatchQueue.main.async { completion(nil) }; return
            }
            let active = jobs.filter { ["queued", "running", "cancelling"].contains($0["status"] as? String ?? "") }
            DispatchQueue.main.async { completion(active) }
        }.resume()
    }

    func applicationShouldTerminate(_ sender: NSApplication) -> NSApplication.TerminateReply {
        if isQuitting { return .terminateNow }
        if terminationPending { return .terminateCancel }
        if !downloads.isEmpty || savePanels > 0 {
            showError("A download is still being saved", "Finish or cancel the download before quitting DynaMol.")
            return .terminateCancel
        }
        terminationPending = true
        // An asynchronous reply also covers Quit from the Dock and system logout.
        DispatchQueue.main.async { [weak self] in self?.confirmQuit() }
        return .terminateLater
    }

    private func confirmQuit() {
        if serviceState() == nil && child?.isRunning == true { finishQuit(keepService: false); return }
        activeJobs { [weak self] jobs in
            guard let self = self else { return }
            if let jobs = jobs, jobs.isEmpty { self.finishQuit(keepService: false); return }
            if jobs == nil && self.serviceState() == nil && self.child == nil { self.finishQuit(keepService: false); return }
            self.openWorkspace()
            let alert = NSAlert()
            alert.messageText = jobs == nil ? "The local service could not confirm job status" : "Simulation work is still running"
            alert.informativeText = jobs == nil
                ? "You can keep DynaMol open or quit while leaving the local service running. Reopen DynaMol to check your jobs."
                : "\(jobs!.count) background job(s) are active. Keep DynaMol open to watch progress, stop the jobs before quitting, or leave them running and reopen DynaMol later. Stopped MD jobs can resume from saved production checkpoints; preparation must be started again."
            alert.addButton(withTitle: "Keep DynaMol Open")
            if jobs != nil { alert.addButton(withTitle: "Stop Jobs and Quit") }
            alert.addButton(withTitle: "Keep Jobs Running and Quit")
            let answer = alert.runModal()
            if answer == .alertFirstButtonReturn { self.cancelQuit() }
            else if jobs != nil && answer == .alertSecondButtonReturn { self.stopJobs(jobs!) }
            else { self.finishQuit(keepService: true) }
        }
    }

    private func cancelQuit() { terminationPending = false; NSApp.reply(toApplicationShouldTerminate: false) }
    private func stopJobs(_ jobs: [[String: Any]]) {
        guard let state = serviceState(), let address = state["url"] as? String, let base = URL(string: address) else { cancelQuit(); return }
        window.subtitle = "Stopping background jobs…"
        let group = DispatchGroup()
        for job in jobs {
            guard let id = job["id"] as? String else { continue }
            var request = URLRequest(url: base.appendingPathComponent("api/jobs").appendingPathComponent(id).appendingPathComponent("cancel"))
            request.httpMethod = "POST"; request.timeoutInterval = 5
            group.enter()
            URLSession.shared.dataTask(with: request) { _, _, _ in group.leave() }.resume()
        }
        group.notify(queue: .main) { [weak self] in self?.waitForStoppedJobs(deadline: Date().addingTimeInterval(45)) }
    }
    private func waitForStoppedJobs(deadline: Date) {
        activeJobs { [weak self] jobs in
            guard let self = self else { return }
            if let jobs = jobs, jobs.isEmpty { self.finishQuit(keepService: false); return }
            if Date() < deadline {
                DispatchQueue.main.asyncAfter(deadline: .now() + 1) { self.waitForStoppedJobs(deadline: deadline) }
            } else {
                self.window.subtitle = ""
                self.cancelQuit()
                self.showError("DynaMol is waiting for the engine to stop", "The workspace remains open so you can inspect job status. Saved files and complete production checkpoints are retained.")
            }
        }
    }

    private func finishQuit(keepService: Bool) {
        isQuitting = true
        leaveServiceRunning = keepService
        startupTimer?.invalidate()
        if keepService { NSApp.reply(toApplicationShouldTerminate: true); return }
        if let process = child, process.isRunning {
            process.terminate()
            DispatchQueue.global().async {
                // Give the local service time to close its files and release its launch lock.
                let deadline = Date().addingTimeInterval(18)
                while process.isRunning && Date() < deadline { Thread.sleep(forTimeInterval: 0.1) }
                DispatchQueue.main.async { NSApp.reply(toApplicationShouldTerminate: true) }
            }
        } else if let object = serviceState(), let address = object["control_url"] as? String,
                  let token = object["control_token"] as? String, let url = URL(string: address),
                  url.scheme == "http", url.host == "127.0.0.1", url.port != nil {
            var request = URLRequest(url: url)
            request.httpMethod = "POST"
            request.setValue(token, forHTTPHeaderField: "X-DynaMol-Key")
            request.timeoutInterval = 3
            URLSession.shared.dataTask(with: request) { _, _, _ in
                DispatchQueue.main.async { NSApp.reply(toApplicationShouldTerminate: true) }
            }.resume()
        } else { NSApp.reply(toApplicationShouldTerminate: true) }
    }
    func applicationWillTerminate(_ notification: Notification) {
        if !leaveServiceRunning, let process = child, process.isRunning { process.terminate() }
    }

    private func localNavigation(_ url: URL) -> Bool {
        if url.absoluteString == "about:blank" { return true }
        if url.scheme == "blob", let nested = URL(string: String(url.absoluteString.dropFirst(5))) { return localNavigation(nested) }
        guard url.scheme == "http", url.host == "127.0.0.1" else { return false }
        return [stateFile("startup.json"), serviceState()].compactMap { $0?["url"] as? String }
            .compactMap(URL.init(string:)).contains { $0.scheme == url.scheme && $0.host == url.host && $0.port == url.port }
    }
    func webView(_ webView: WKWebView, decidePolicyFor navigationAction: WKNavigationAction, decisionHandler: @escaping (WKNavigationActionPolicy) -> Void) {
        guard let url = navigationAction.request.url else { decisionHandler(.cancel); return }
        guard localNavigation(url) else {
            if url.scheme == "https", navigationAction.navigationType == .linkActivated { NSWorkspace.shared.open(url) }
            decisionHandler(.cancel)
            return
        }
        decisionHandler(navigationAction.shouldPerformDownload ? .download : .allow)
    }
    func webView(_ webView: WKWebView, decidePolicyFor navigationResponse: WKNavigationResponse, decisionHandler: @escaping (WKNavigationResponsePolicy) -> Void) {
        let attachment = (navigationResponse.response as? HTTPURLResponse)?.value(forHTTPHeaderField: "Content-Disposition")?.lowercased().hasPrefix("attachment") == true
        decisionHandler(attachment || !navigationResponse.canShowMIMEType ? .download : .allow)
    }
    func webView(_ webView: WKWebView, createWebViewWith configuration: WKWebViewConfiguration, for navigationAction: WKNavigationAction, windowFeatures: WKWindowFeatures) -> WKWebView? {
        if let url = navigationAction.request.url {
            if localNavigation(url) { webView.load(navigationAction.request) }
            else if url.scheme == "https" { NSWorkspace.shared.open(url) }
        }
        return nil
    }
    func webView(_ webView: WKWebView, didFailProvisionalNavigation navigation: WKNavigation!, withError error: Error) {
        if (error as NSError).code != NSURLErrorCancelled && !isQuitting {
            showError("The workspace could not load", "\(error.localizedDescription) Use View → Reload Workspace to try again.")
        }
    }
    func webViewWebContentProcessDidTerminate(_ webView: WKWebView) {
        showError("The molecular viewer needs to reload", "Saved projects and simulation jobs are safe. The viewer will now reload; any unsaved view changes may need to be made again.")
        webView.reload()
    }
    func webView(_ webView: WKWebView, runOpenPanelWith parameters: WKOpenPanelParameters, initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping ([URL]?) -> Void) {
        let panel = NSOpenPanel()
        panel.title = "Open molecular files"
        panel.canChooseFiles = true
        panel.canChooseDirectories = parameters.allowsDirectories
        panel.allowsMultipleSelection = parameters.allowsMultipleSelection
        panel.beginSheetModal(for: window) { answer in completionHandler(answer == .OK ? panel.urls : nil) }
    }
    func webView(_ webView: WKWebView, runJavaScriptAlertPanelWithMessage message: String, initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping () -> Void) {
        let alert = NSAlert(); alert.messageText = "DynaMol"; alert.informativeText = message
        alert.beginSheetModal(for: window) { _ in completionHandler() }
    }
    func webView(_ webView: WKWebView, runJavaScriptConfirmPanelWithMessage message: String, initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping (Bool) -> Void) {
        let alert = NSAlert(); alert.messageText = "DynaMol"; alert.informativeText = message
        alert.addButton(withTitle: "OK"); alert.addButton(withTitle: "Cancel")
        alert.beginSheetModal(for: window) { result in completionHandler(result == .alertFirstButtonReturn) }
    }

    func webView(_ webView: WKWebView, navigationAction: WKNavigationAction, didBecome download: WKDownload) { retain(download) }
    func webView(_ webView: WKWebView, navigationResponse: WKNavigationResponse, didBecome download: WKDownload) { retain(download) }
    private func retain(_ download: WKDownload) { downloads[ObjectIdentifier(download)] = download; download.delegate = self }
    func download(_ download: WKDownload, decideDestinationUsing response: URLResponse, suggestedFilename: String, completionHandler: @escaping (URL?) -> Void) {
        // WebKit needs a new destination immediately. Stage first so revoked Blob URLs
        // and an existing filename in the eventual Save dialog remain safe.
        let folder = FileManager.default.temporaryDirectory.appendingPathComponent("DynaMol-download-" + UUID().uuidString, isDirectory: true)
        do {
            try FileManager.default.createDirectory(at: folder, withIntermediateDirectories: true)
            let basename = URL(fileURLWithPath: suggestedFilename).lastPathComponent
            let file = folder.appendingPathComponent(basename.isEmpty || basename == "." || basename == "/" ? "DynaMol-download" : basename)
            downloadFiles[ObjectIdentifier(download)] = file
            window.subtitle = "Downloading \(file.lastPathComponent)…"
            completionHandler(file)
        } catch { completionHandler(nil); showError("The download could not start", error.localizedDescription) }
    }
    func downloadDidFinish(_ download: WKDownload) {
        let key = ObjectIdentifier(download)
        downloads.removeValue(forKey: key)
        guard let temporary = downloadFiles.removeValue(forKey: key) else { return }
        window.subtitle = ""
        savePanels += 1
        let panel = NSSavePanel()
        panel.title = "Save DynaMol export"
        panel.nameFieldStringValue = temporary.lastPathComponent
        panel.canCreateDirectories = true
        panel.beginSheetModal(for: window) { [weak self] answer in
            guard let self = self else { return }
            self.savePanels -= 1
            guard answer == .OK, let destination = panel.url else {
                try? FileManager.default.removeItem(at: temporary.deletingLastPathComponent())
                return
            }
            do {
                if FileManager.default.fileExists(atPath: destination.path) {
                    _ = try FileManager.default.replaceItemAt(destination, withItemAt: temporary)
                } else { try FileManager.default.moveItem(at: temporary, to: destination) }
                try? FileManager.default.removeItem(at: temporary.deletingLastPathComponent())
            } catch {
                // Keep the completed file reachable when a destination volume is unavailable.
                let recovery = self.support.appendingPathComponent("Recovered Downloads", isDirectory: true)
                try? FileManager.default.createDirectory(at: recovery, withIntermediateDirectories: true)
                let saved = recovery.appendingPathComponent(UUID().uuidString + "-" + temporary.lastPathComponent)
                do {
                    try FileManager.default.moveItem(at: temporary, to: saved)
                    try? FileManager.default.removeItem(at: temporary.deletingLastPathComponent())
                    self.showError("The export could not be saved there", "\(error.localizedDescription)\nA recovery copy is available at \(saved.path).")
                } catch {
                    self.showError("The export could not be saved there", "\(error.localizedDescription)\nThe completed download is still available at \(temporary.path).")
                }
            }
        }
    }
    func download(_ download: WKDownload, didFailWithError error: Error, resumeData: Data?) {
        let key = ObjectIdentifier(download)
        downloads.removeValue(forKey: key)
        if let file = downloadFiles.removeValue(forKey: key) { try? FileManager.default.removeItem(at: file.deletingLastPathComponent()) }
        window.subtitle = ""
        if (error as NSError).code != NSURLErrorCancelled { showError("The download did not finish", error.localizedDescription) }
    }
    private func showError(_ title: String, _ detail: String) {
        let alert = NSAlert(); alert.messageText = title; alert.informativeText = detail
        if let window = window { window.makeKeyAndOrderFront(nil); alert.beginSheetModal(for: window) }
        else { alert.runModal() }
    }
}

let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.setActivationPolicy(.regular)
app.run()
