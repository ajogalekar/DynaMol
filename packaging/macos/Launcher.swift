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
    @objc func quit() {
        let alert = NSAlert()
        alert.messageText = "Quit DynaMol?"
        alert.informativeText = "The local viewer service will close. Simulation workers already running keep their saved progress and can be reopened the next time you launch DynaMol."
        alert.addButton(withTitle: "Quit")
        alert.addButton(withTitle: "Keep open")
        if alert.runModal() == .alertFirstButtonReturn {
            if let process = child, process.isRunning {
                process.terminate()
                NSApplication.shared.terminate(nil)
            } else {
                // Reopened launchers reuse the service through a per-run local
                // control token, never by signalling a potentially stale PID.
                let state = support.appendingPathComponent("service.json")
                if let data = try? Data(contentsOf: state),
                   let object = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
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
        }
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
