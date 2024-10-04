use eframe::egui;
use lute::prelude::*;

#[cfg(not(target_family = "wasm"))]
fn main() -> eframe::Result {
    let mut arena = Arena::<u32>::new();
    let mut current = 0;
    let mut smiles_input = String::new();
    let mut smiles_popup_focus = true;
    eframe::run_simple_native(
        "Lute GUI",
        eframe::NativeOptions::default(),
        move |ctx, _frame| {
            egui::SidePanel::left("sidebar").show(ctx, |ui| {
                ui.horizontal_top(|ui| {
                    ui.label("Fragments");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                        if ui.button("C").clicked() {
                            arena = Arena::new();
                        }
                        {
                            let popup_id = ui.make_persistent_id("load_smiles");
                            let response = ui.button("+");
                            if response.clicked() {
                                ui.memory_mut(|mem| mem.open_popup(popup_id));
                                smiles_popup_focus = true;
                            }
                            egui::popup_above_or_below_widget(
                                ui,
                                popup_id,
                                &response,
                                egui::AboveOrBelow::Below,
                                egui::PopupCloseBehavior::IgnoreClicks,
                                |ui| {
                                    ui.set_min_width(100.0);
                                    let editor = ui.add(
                                        egui::TextEdit::singleline(&mut smiles_input)
                                            .font(egui::TextStyle::Monospace)
                                            .return_key(None),
                                    );
                                    if smiles_popup_focus {
                                        editor.request_focus();
                                        smiles_popup_focus = false;
                                    }
                                    ui.with_layout(
                                        egui::Layout::top_down_justified(egui::Align::Center),
                                        |ui| {
                                            let load = ui.button("Load");
                                            let cancel = ui.button("Cancel");
                                            if load.clicked()
                                                || (editor.has_focus()
                                                    && ctx
                                                        .input(|i| i.key_pressed(egui::Key::Enter)))
                                            {
                                                let graph = SmilesParser::new(&smiles_input)
                                                    .parse()
                                                    .expect("TODO");
                                                let idx = arena.insert_mol(&graph);
                                                current = idx.0;
                                                smiles_input.clear();
                                                ui.memory_mut(|mem| mem.close_popup());
                                            } else if cancel.clicked()
                                                || (editor.has_focus()
                                                    && ctx.input(|i| {
                                                        i.key_pressed(egui::Key::Escape)
                                                    }))
                                            {
                                                smiles_input.clear();
                                                ui.memory_mut(|mem| mem.close_popup());
                                            }
                                        },
                                    );
                                },
                            );
                        }
                    });
                });
                ui.separator();
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.with_layout(egui::Layout::top_down_justified(egui::Align::Min), |ui| {
                        for n in 0..arena.expose_frags().len() {
                            let n = n as u32;
                            let button = ui.button(format!(
                                "{n}: {}",
                                arena.molecule(n.into()).smiles(Default::default())
                            ));
                            if button.clicked() {
                                current = n;
                            }
                        }
                    });
                });
            });
            egui::CentralPanel::default().show(ctx, |ui| {
                if !arena.expose_frags().is_empty() {
                    ui.centered_and_justified(|ui| {
                        let mol = arena.molecule(current.into());

                        let data = SvgFormatter {
                            graph: mol,
                            mode: lute::disp::svg::FormatMode::LegacyH,
                        }
                        .render(None);
                        let img = egui::ColorImage {
                            size: [data.width() as _, data.height() as _],
                            pixels: bytemuck::allocation::cast_vec::<u8, egui::Color32>(
                                data.take(),
                            ),
                        };
                        let texture = ctx.load_texture(
                            format!("mol:{}", mol.lightweight_smiles()),
                            img,
                            Default::default(),
                        );
                        ui.add(egui::Image::new((texture.id(), texture.size_vec2())));
                        // TODO: make this cache the image
                    });

                    let mol = arena.molecule(current.into());

                    ui.with_layout(egui::Layout::bottom_up(egui::Align::Max), |ui| {
                        ui.label(format!(
                            "Fast SMILES: {}",
                            mol.smiles(SmilesConfig::fast_roundtrip()),
                        ));
                        ui.label(format!("Canon SMILES: {}", mol.smiles(SmilesConfig::new())));
                        ui.label(format!("Atomic mass: {:.3}u", mol.mass()));
                    });
                }
            });
        },
    )
}
