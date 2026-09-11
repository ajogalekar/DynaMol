import AppKit
import Foundation

final class AppDelegate: NSObject, NSApplicationDelegate {
    private var child: Process?
    private var item: NSStatusItem?
    private let support: URL = {
        let arguments = CommandLine.arguments
        if let index = arguments.firstIndex(of: "--home"), index + 1 < arguments.count {
            return URL(fileURLWithPath: arguments[index + 1], isDirectory: true)
        }
        return FileManager.default.homeDirectoryForCurrentUser.appendingPathComponent("Library/Application Support/DynaMol")
    }()
    func applicationDidFinishLaunching(_ notification: Notification) {
        item = NSStatusBar.system.statusItem(withLength: NSStatusItem.variableLength)
        item?.button?.title = "DynaMol"
        let menu = NSMenu()
        menu.addItem(withTitle: "Open DynaMol", action: #selector(openWorkspace), keyEquivalent: "o").target = self
        menu.addItem(withTitle: "Show workspace files", action: #selector(showFiles), keyEquivalent: "").target = self
        menu.addItem(NSMenuItem.separator())
        menu.addItem(withTitle: "Quit DynaMol", action: #selector(quit), keyEquivalent: "q").target = self
        item?.menu = menu
        start()
    }
    private func start() {
        guard let resources = Bundle.main.resourceURL else { return }
        let process = Process()
        process.executableURL = resources.appendingPathComponent("python/bin/python3.12")
        process.arguments = ["-I", "-B", resources.appendingPathComponent("bootstrap.py").path] + Array(CommandLine.arguments.dropFirst())
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
                if process.terminationStatus != 0 {
                    let alert = NSAlert()
                    alert.messageText = "DynaMol could not start"
                    alert.informativeText = "Your molecules are safe. See Library/Application Support/DynaMol/Logs for details."
                    alert.runModal()
                }
                self?.child = nil
            }
        }
        do { try process.run(); child = process }
        catch {
            let alert = NSAlert(); alert.messageText = "DynaMol could not start"; alert.informativeText = error.localizedDescription; alert.runModal()
        }
    }
    @objc func openWorkspace() {
        let state = support.appendingPathComponent("service.json")
        if let data = try? Data(contentsOf: state), let object = try? JSONSerialization.jsonObject(with: data) as? [String: Any], let address = object["url"] as? String, let url = URL(string: address) {
            NSWorkspace.shared.open(url)
        } else if let data = try? Data(contentsOf: support.appendingPathComponent("startup.json")),
                  let object = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
                  let address = object["url"] as? String,
                  let url = URL(string: address), url.host == "127.0.0.1" {
            NSWorkspace.shared.open(url)
        } else if child == nil { start() }
    }
    @objc func showFiles() {
        NSWorkspace.shared.open(support.appendingPathComponent("Workspace"))
    }
    private func serviceState() -> [String: Any]? {
        guard let data = try? Data(contentsOf: support.appendingPathComponent("service.json")),
              let object = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
              let address = object["url"] as? String, let url = URL(string: address), url.host == "127.0.0.1" else { return nil }
        return object
    }
    private func activeJobs(_ completion: @escaping ([[String: Any]]?) -> Void) {
        guard let state = serviceState(), let address = state["url"] as? String,
              let url = URL(string: address + "api/jobs") else { completion(nil); return }
        var request = URLRequest(url: url); request.timeoutInterval = 3
        URLSession.shared.dataTask(with: request) { data, response, _ in
            guard (response as? HTTPURLResponse)?.statusCode == 200, let data = data,
                  let jobs = try? JSONSerialization.jsonObject(with: data) as? [[String: Any]] else {
                DispatchQueue.main.async { completion(nil) }; return
            }
            let active = jobs.filter { ["queued", "running", "cancelling"].contains($0["status"] as? String ?? "") }
            DispatchQueue.main.async { completion(active) }
        }.resume()
    }
    @objc func quit() {
        activeJobs { [weak self] jobs in
            guard let self = self else { return }
            if let jobs = jobs, jobs.isEmpty { self.shutdown(); return }
            let alert = NSAlert()
            alert.messageText = jobs == nil ? "Close the DynaMol viewer?" : "Simulation work is still running"
            alert.informativeText = jobs == nil
                ? "The local service could not confirm job status. Closing the viewer leaves any background workers running. Reopen DynaMol to check their progress."
                : "\(jobs!.count) background job(s) are active. You can keep them running after the viewer closes, or stop them first. MD jobs can resume from their last saved production checkpoint; preparation work must be started again."
            alert.addButton(withTitle: "Quit viewer, keep running")
            if jobs != nil { alert.addButton(withTitle: "Stop jobs and quit") }
            alert.addButton(withTitle: "Keep DynaMol open")
            let answer = alert.runModal()
            if answer == .alertFirstButtonReturn { self.shutdown() }
            else if jobs != nil && answer == .alertSecondButtonReturn { self.stopJobs(jobs!) }
        }
    }
    private func stopJobs(_ jobs: [[String: Any]]) {
        guard let state = serviceState(), let address = state["url"] as? String else { return }
        item?.button?.title = "DynaMol · stopping…"
        let group = DispatchGroup()
        for job in jobs {
            guard let id = job["id"] as? String, let url = URL(string: address + "api/jobs/" + id + "/cancel") else { continue }
            var request = URLRequest(url: url); request.httpMethod = "POST"; request.timeoutInterval = 5
            group.enter()
            URLSession.shared.dataTask(with: request) { _, _, _ in group.leave() }.resume()
        }
        group.notify(queue: .main) { [weak self] in self?.waitForStoppedJobs(deadline: Date().addingTimeInterval(45)) }
    }
    private func waitForStoppedJobs(deadline: Date) {
        activeJobs { [weak self] jobs in
            guard let self = self else { return }
            if let jobs = jobs, jobs.isEmpty { self.shutdown(); return }
            if Date() < deadline {
                DispatchQueue.main.asyncAfter(deadline: .now() + 1) { self.waitForStoppedJobs(deadline: deadline) }
                return
            }
            self.item?.button?.title = "DynaMol"
            let alert = NSAlert(); alert.messageText = "DynaMol is waiting for the engine to stop"
            alert.informativeText = "The viewer remains open so you can inspect the job status. Saved files and complete production checkpoints are retained."
            alert.addButton(withTitle: "Open DynaMol"); alert.runModal(); self.openWorkspace()
        }
    }
    private func shutdown() {
        if let process = child, process.isRunning {
            process.terminate()
            NSApplication.shared.terminate(nil)
        } else if let object = serviceState(),
                  let address = object["control_url"] as? String,
                  let token = object["control_token"] as? String,
                  let url = URL(string: address), url.host == "127.0.0.1" {
            var request = URLRequest(url: url)
            request.httpMethod = "POST"
            request.setValue(token, forHTTPHeaderField: "X-DynaMol-Key")
            request.timeoutInterval = 2
            URLSession.shared.dataTask(with: request) { _, _, _ in
                DispatchQueue.main.async { NSApplication.shared.terminate(nil) }
            }.resume()
        } else { NSApplication.shared.terminate(nil) }
    }
    func applicationWillTerminate(_ notification: Notification) {
        if let process = child, process.isRunning { process.terminate() }
    }
    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { false }
}

let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.setActivationPolicy(.accessory)
app.run()
